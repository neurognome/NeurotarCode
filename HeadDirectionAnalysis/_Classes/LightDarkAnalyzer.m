classdef LightDarkAnalyzer < handle
    properties
        light_hda
        dark_hda
        
        light_data
        dark_data
        
        clust_id
        n_clusters
    end
    
    methods
        function obj = LightDarkAnalyzer(light_hda, dark_hda)
            obj.light_hda = light_hda; % Here we inject the two separate HeadDirectionAnalyzers into this class for access
            obj.dark_hda = dark_hda;
            
            obj.light_data = obj.light_hda.getBinnedData();
            obj.dark_data = obj.dark_hda.getBinnedData();
            
            % REmove junk cells first
            
            %
            %                     obj.light_data = obj.light_data(obj.light_hda.is_head_direction, :);
            %             obj.dark_data = obj.dark_data(obj.light_hda.is_head_direction, :);
        end
        
        function classifyResponses(obj, check_flag, light_tuning, dark_tuning)
            if nargin < 2 || isempty(check_flag)
                check_flag = true;
            end

            if nargin < 3 || isempty(light_tuning)
                light_tuning = obj.light_data;
            end
            
            if nargin < 4 || isempty(dark_tuning)
                dark_tuning = obj.dark_data;
            end
            
            combined_tuning = mean(cat(3, light_tuning, dark_tuning), 3);
            data = zeros(size(light_tuning, 1), (size(dark_tuning, 2) + size(light_tuning, 2)));
            
            % Preprocess and fit the data 
            parfor c = 1:size(light_tuning, 1) % parallel for loop wows peed
                % Find the proper offset
                fprintf('Fitting gaussian for cell %d/%d \n', c, size(light_tuning, 1))
                bound = @(x, bl, bu) min(max(x, bl), bu);
                modelfun = @(b, x) b(1) + ...
                max([0, b(2)]) * exp(-(x(:, 1) - bound(b(3), 0, x(end, 1))) .^ 2 ./ (2 * b(4) .^ 2)) + ...
                max([0, b(5)]) * exp(-(x(:, 1) - bound(b(6), 0, x(end, 1))) .^ 2 ./ (2 * b(4) .^ 2)); % Based on Michael natneuro paper

                coeffs = obj.fitDoubleGaussian(combined_tuning(c, :), modelfun);

                [coeffs_light, rsqr] = obj.fitDoubleGaussian(light_tuning(c, :), modelfun);

                % Baseline shift
                if rsqr > 0
                    light = (circshift(light_tuning(c, :), size(light_tuning, 2)/4 - round(coeffs(3)))) - coeffs_light(1); % baseline sub
                else
                    light = (circshift(light_tuning(c, :), size(light_tuning, 2)/4 - round(coeffs(3)))) - mean(light_tuning(c, :));
                end
                
                [coeffs_dark, rsqr] = obj.fitDoubleGaussian(dark_tuning(c, :), modelfun);
                if rsqr > 0
                    dark = (circshift(dark_tuning(c, :), size(light_tuning, 2)/4 - round(coeffs(3)))) - coeffs_dark(1);
                else
                    dark = (circshift(dark_tuning(c, :), size(dark_tuning, 2)/4 - round(coeffs(3)))) - mean(dark_tuning(c, :));
                end

                % Put the processed data together for clustering
                data(c, :) = cat(2, rescale(light), rescale(dark));
            end
            % data = data(obj.light_hda.is_head_direction, :);

            % Run clustering
            max_clust_flag = false;
            while true
                obj.clust_id = obj.cluster(data, true, max_clust_flag); % PCA clustering
                obj.n_clusters = max(unique(obj.clust_id));
                
                if check_flag
                    obj.checkClusteringPerformance(data)
                end
                happy = input('Happy? (0 no 1 yes): ');
                if happy
                    break
                else
                    max_clust_flag = true;
                end
            end
        end
        
        function viewClusteredResponses(obj)
            unshifted_data = zeros(size(obj.light_data, 1), 2 * size(obj.light_data, 2));
            for c = 1:size(obj.light_data, 1)
                light = zscore(obj.light_data(c, :));
                
                dark = zscore(obj.dark_data(c, :));
                unshifted_data(c, :) = (cat(2, light,  dark));
            end
            
            grouped_data = cell(1, max(unique(obj.clust_id)));
            ct = 1;
            for c = unique(obj.clust_id)'
                curr = unshifted_data(obj.clust_id == c, :);
                [~, max_idx] = max(curr(:, 1:size(unshifted_data, 2) / 2), [], 2);
                [~, sort_idx] = sort(max_idx);
                grouped_data{ct} = curr(sort_idx, :);
                ct = ct + 1;
            end
            
            figure
            for i_group = 1:length(grouped_data)
                imagesc(grouped_data{i_group});
                title(sprintf('Group #%d', i_group))
                colormap(flipud(bone))
                pause
            end
        end
        
        function out = getClusters(obj)
            out = obj.clust_id;
        end
        
        function clustID = cluster(obj, data, pca_flag, max_clust_flag)
            % From Michael 20Dec2019
            if nargin < 3 || isempty(pca_flag)
                pca_flag = true;
            end

            varExpl = 0.98;    % criteria for variance explained (deteremines # of PCs)
            
            numSamp = size(data, 2);
            numCells = size(data, 1);
            
            %% PCA
            if pca_flag
                [~,score,latent] = pca(data);
                PC_var = cumsum(latent/sum(latent));
                numPCs = find(PC_var > varExpl,1);
                PCA_resp = score(:, 1:numPCs);
            else
                PCA_resp = data;
            end

            %% Cluster responses

            max_clust_flag = false;
            Z = linkage(PCA_resp, 'ward', 'euclidean');
            f = figure;
            if max_clust_flag
                [~, ~, cluster_order] = dendrogram(Z, numCells)
                set(gcf,'color',[1 1 1])
                numClust = input('How many clusters? ');    % Determine cluster number on dendrogram
                clustID = cluster(Z,'maxclust',numClust);
            else
                [~, ~, cluster_order] = dendrogram(Z, numCells, 'ColorThreshold', 9)
                set(gcf,'color',[1 1 1])
                clustID = cluster(Z, 'cutoff', 9, 'criterion', 'distance');
                numClust = numel(unique(clustID)); 
            end
            order = unique(clustID(cluster_order), 'stable'); % don't sort
            close(f)
            
            % Display the cluster ID
            f = figure;
            dendrogram(Z, numCells, 'ColorThreshold', 9)
            hold on
            line([xlim], [9, 9]);
            labels = str2num(xticklabels);
            ct = 1;
            for ii = order'
                mid_val = round(median(find(clustID(cluster_order) == ii)));
                for c = 1:numCells
                    if labels(c) == mid_val
                        disp(labels)
                        tick_pos(ct) = labels(c);
                        ct = ct + 1;
                    end
                end
            end
            xticks(tick_pos)
            xticklabels(order)
            xlabel('Cluster ID')

            %% Recluster by correlation to mean response (to clean up sorting)
            clusterMeanTraces = zeros(numClust, numSamp);
            for c = 1:numClust % calculate average traces
                clusterMeanTraces(c,:) = mean(data(clustID == c,:));
            end

            for n = 1:numCells % assign ased on CC to average traces
                curr_trace = data(n, :);
                for c = 1:numClust
                    CC = corrcoef(curr_trace, clusterMeanTraces(c,:));
                    clustCC(c) = CC(2);
                end
                [~, clustID(n)] = max(clustCC);
            end


        end
        
        function checkClusteringPerformance(obj, data)
            obj.returnClusteringMetrics()
            obj.plotResults(data);
        end
        
        function returnClusteringMetrics(obj)
            %% Analysis Info
            cells_per_clust = zeros(1, obj.n_clusters);
            ct = 1;
            for c = unique(obj.clust_id)' % recalculate average traces
                cells_per_clust(ct) = sum(obj.clust_id == c);
                ct = ct + 1;
            end
            
            disp(['Cells per cluster: ' num2str(cells_per_clust)])
            
        end
        
        function plotResults(obj, data)
            cluster_mean_traces = zeros(obj.n_clusters, size(data, 2));
            cluster_se = zeros(obj.n_clusters, size(data, 2));
            ct = 1;
            for c = unique(obj.clust_id)'
                cluster_mean_traces(ct, :) = mean(data(obj.clust_id == c, :));
                cluster_se(ct, :) = std(data(obj.clust_id == c, :), 0, 1) / sqrt(sum(obj.clust_id == c));
                ct = ct + 1;
            end
            
            figure
            for c = 1:obj.n_clusters
                legend_text{c} = ['Cluster #' num2str(c)];
            end
            timeVec = [1:size(data, 2)];
            minPlot = floor(min(min(cluster_mean_traces))/5);
            maxPlot = ceil(max(max(cluster_mean_traces))/5);
            hold on
            plot(timeVec,zeros(1,size(data, 2)),'k:','linewidth',2,'HandleVisibility','off')
            plot(zeros(1,size(data, 2)),linspace(minPlot,maxPlot,size(data, 2)),'k:','linewidth',2,'HandleVisibility','off')
            %plot(3*ones(1,numSamp),linspace(minPlot,maxPlot,numSamp),'k:','linewidth',2,'HandleVisibility','off')
            %plot(6*ones(1,numSamp),linspace(minPlot,maxPlot,numSamp),'k:','linewidth',2,'HandleVisibility','off')
            
            plot(timeVec,cluster_mean_traces','linewidth',3)
            legend(legend_text)
            xlabel('Time (sec)')
            ylabel('DF/F')
            set(gcf,'color',[1 1 1])
            saveas(gcf,'Cluster_AverageTraces','fig')
            
            subplotIndices = [1 2 7 8];
            for c = 1:obj.n_clusters
                if rem(c-1,4)==0
                    figure
                end
                currMatrix = data(obj.clust_id==c,:);
                currMatrix = currMatrix-repmat(min(currMatrix,[],2),[1 size(data, 2)]);
                currMatrix = currMatrix./repmat(max(currMatrix,[],2),[1 size(data, 2)]);
                currMatrix = currMatrix(randperm(size(currMatrix,1)),:);
                cluster_SE_lower = cluster_mean_traces(c,:)-cluster_se(c,:);
                cluster_SE_upper = cluster_mean_traces(c,:)+cluster_se(c,:);
                currSubplotIdx = subplotIndices(rem(c-1,4)+1);
                
                subplot(6,2,[currSubplotIdx currSubplotIdx+2])
                imagesc(timeVec, [1:size(currMatrix,1)],currMatrix)
                title(['Cluster #' num2str(c) ', n = ' num2str(size(currMatrix,1))])
                
                subplot(6, 2, [currSubplotIdx+4])
                hold on
                % plot(timeVec(1:numSamp/2),zeros(1,numSamp/2),'k:','linewidth',2)
                %  plot(zeros(1,numSamp),linspace(minPlot,maxPlot,numSamp),'k:','linewidth',2)
                timeVecShort = timeVec(1:size(data, 2)/2);
                
                se_low = cluster_SE_lower(1:size(data, 2)/2);
                se_hi = cluster_SE_upper(1:size(data, 2)/2);
                fill([timeVecShort, fliplr(timeVecShort)], [se_low, fliplr(se_hi)], [0 0 1],'facealpha', 0.5,'linestyle','none')
                se_low = cluster_SE_lower(size(data, 2)/2 + 1:end);
                se_hi = cluster_SE_upper(size(data, 2)/2 + 1:end);
                fill([timeVecShort, fliplr(timeVecShort)], [se_low, fliplr(se_hi)], [1 0 0],'facealpha', 0.5,'linestyle','none')
                xlabel('Time (sec)')
                ylabel('DF/F')
                if rem(c-1,4)==3 || c==obj.n_clusters
                    set(gcf,'color',[1 1 1])
                    set(gcf,'position',[850 80 550 900])
                    figCount = ceil(c/4);
                    %saveas(gcf,['Cluster_Responses_' num2str(figCount)],'fig')
                end
                
            end
        end
        
        function [coefficients, fit_quality] = fitDoubleGaussian(obj, data, model)
            % The plan is to use a mixture of gaussian models to get a better idea of the bilobality stuff
            % Prepare options for linear fitting
            opts = statset('nlinfit');
            opts.FunValCheck = 'on';
            %opts.Display = 'final';
            coefficients = zeros(size(data, 1), 6);
            fit_quality = zeros(size(data, 1), 1);
            for c = 1:size(data, 1)
                curr_data = data(c, :);
                
                beta0 = obj.prepareBeta0(curr_data);
                tbl = table([1:size(data, 2)]', data(c, :)');
                try
                    mdl = fitnlm(tbl, model, beta0, 'Options', opts);
                    coefficients(c, :) = mdl.Coefficients{:, 'Estimate'};
                    fit_quality(c) = mdl.Rsquared.Adjusted;
                catch
                    coefficients(c, :) = NaN;
                    fit_quality(c) = 0;
                end
                
            end
        end
        
        function beta0 = prepareBeta0(obj, data)
            [sorted_values, sorted_idx] = sort(data, 'descend');
            mag1 = max([sorted_values(1), 1]);
            loc1 = sorted_idx(1);
            std1 = 2;
            
            % find the second peak
            wrapN = @(x, N) (1 + mod(x-1, N));
            mag2 = mag1 / 2; % Random guess
            loc2 = wrapN(loc1 + length(data) / 2, length(data));
            
            beta0 = [0, mag1, loc1, std1, mag2, loc2];
        end

        function correctDarkDrift(obj)

            %% this is filthy, super spaghetti... but it works.. for now lol
            n_segments = 5;
            
            data = obj.dark_hda.data.spikes;
            heading = obj.dark_hda.floating.heading;
            segmented_data = obj.segmentData(data, n_segments);
            segmented_heading = obj.segmentData(heading', n_segments);
            
            tuning_curve = zeros(size(obj.dark_data, 1), size(obj.dark_data, 2), n_segments);
            for i_segment = 1:n_segments
                tuning_curve(:, :, i_segment) = obj.binData([], segmented_data{i_segment}, segmented_heading{i_segment});
            end
            bound = @(x, bl, bu) min(max(x, bl), bu);
            modelfun = @(b, x) b(1) + ...
            max([0, b(2)]) * exp(-(x(:, 1) - bound(b(3), 0, x(end, 1))) .^ 2 / (2 * b(4) .^ 2)) + ...
                    max([0, b(5)]) * exp(-(x(:, 1) - bound(b(6), 0, x(end, 1))) .^ 2 / (2 * b(4) .^ 2)); % Based on Michael natneuro paper

            % Fitting peak for each sigment
            disp('Fitting peaks for all cells across each segment...')
            for i_cell = 1:size(tuning_curve, 1)
                for i_segment = 1:n_segments
                    coeffs = obj.fitDoubleGaussian(tuning_curve(i_cell, :, i_segment), modelfun);
                    peak(i_cell, i_segment) = coeffs(:, 3);
                end
            end
            
            % remove max and min peaks, assume no good
            peak(peak == 1 | peak == 60) = NaN;
            light_peak = obj.light_hda.direction_info(:, 1) .* 60 / (2*pi);
            
            drift = light_peak - peak;
            
            % across the population drift
            pop_drift = nanmedian(drift);
            
            shifted_tuning_curves = zeros(size(tuning_curve));
            for i_segment = 1:n_segments
                shifted_tuning_curves(:, :, i_segment) = circshift(tuning_curve(:, :, i_segment),...
                    round(pop_drift(i_segment)), 2);
            end
            
            obj.dark_data = mean(shifted_tuning_curves, 3);
            
            
            
            % Adjust each to the mean?
%             mean_peak = round(mean(peak, 2));
%             rounded_peaks = round(peak);
%             
%             for i_cell = 1:size(tuning_curve, 1)
%                 for i_seg = 1:n_segments
%                       curr_trace = squeeze(tuning_curve(i_cell, :, i_seg));
%                       temp(:, i_seg) = circshift(curr_trace, mean_peak(i_cell) - rounded_peaks(i_cell, i_seg));
%                 end
%                 
%                 obj.dark_data(i_cell, :) = mean(temp, 2);
%             end

%             for i_cell = 1:size(tuning_curve, 1)
%                 %fit
%                 coeffs = polyfit(1:n_segments, peak(i_cell, :), 1);
%                 slope = coeffs(1);     
%                     for i_seg = 1:n_segments
%                         curr_trace = squeeze(tuning_curve(i_cell, :, i_seg));
%                         temp(:, i_seg) = circshift(curr_trace, round(-(slope * i_seg)));
%                     end
%                     
%                 obj.dark_data(i_cell, :) = mean(temp, 2);
%             end

end


function segmented_data = segmentData(obj, data, n_sections)
    total_dur = size(data, 2);
    section_length = round(total_dur./n_sections);

    segmented_data = cell(n_sections, 1);
            for q = 1:n_sections % Separate timeseries into sections
                segmented_data{q}    = data(:, ...
                    (q - 1) * section_length + 1 : min(q*section_length, size(data, 2)));
            end
        end
        
        function out = binData(obj, bin_width, data, heading, fold_flag)
            % Bins data based on a determined bindwidth
            
            %% added a bunch of stuff in case we want folding.. this is b/c we only want to fold sometimes
            if nargin < 2 || isempty(bin_width)
                bin_width = obj.light_hda.bin_width;
            end
            
            if nargin < 3 ||  isempty(data)
                data = obj.initializeData(obj.data, 'spikes');
            end
            if nargin < 4 || isempty(heading)
                heading = obj.initializeData(obj.floating, 'heading');
            end
            
            if nargin < 5 || isempty(fold_flag)
                fold_flag = false;
            end
            
            % Changed this to be bin_width of 3, and a smoothing filter, a la Giocomo et al 2014, Curr Bio
            if fold_flag
            bin_edges = -360:bin_width:360; % Because alpha is [-180,180];
        else 
            bin_edges = -180:bin_width:180;
        end
            % This is the doubling procedure that's documented in the Jeffrey paper
            %heading = wrapTo180(heading * 2);    
            if fold_flag
            groups = discretize(heading * 2, bin_edges); % doubled
        else
            groups = discretize(heading, bin_edges);
        end
            u_groups = 1:length(bin_edges) - 1; % Changed because sometimes not all groups represented?
            out = zeros(size(data, 1), length(u_groups));
            
            for bin = 1:length(u_groups)
                for c = 1:size(data, 1)
                    out(c, bin) = mean(data(c, groups == u_groups(bin)));
                end
            end
            
            out(isnan(out)) = 0; % NanCheck
            
            % wrap to 2pi
            if fold_flag
                out = mean(cat(3, out(:, 1:120), out(:, 121:end)), 3);
            end

            out = movmean(out, 10, 2); % 10 bins, 5 on each side, = 15 degree on each side, same as Giocomo et al 2014
        end
    end
end