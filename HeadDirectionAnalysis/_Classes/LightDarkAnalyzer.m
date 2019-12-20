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
        end
        
        function classifyResponses(obj, check_flag)
            if nargin < 2 || isempty(check_flag)
                check_flag = true;
            end
            
            data = zeros(size(obj.light_data, 1), 2 * size(obj.light_data, 2));
            for c = 1:size(obj.light_data, 1)
                [~, shift] = max(obj.light_data(c, :));
                light = rescale(circshift(obj.light_data(c, :), 30 - shift));
                
                %  [~, shift] = max(data_d(c, :));
                
                dark = rescale(circshift(obj.dark_data(c, :), 30 - shift));
                data(c, :) = (cat(2, light,  dark));
            end
            
            obj.clust_id = obj.cluster(data);
            obj.n_clusters = max(unique(obj.clust_id));
            
            if check_flag
                obj.checkClusteringPerformance(data)
            end
        end
        
        function viewClusteredResponses(obj)
            unshifted_data = zeros(size(obj.light_data, 1), 2 * size(obj.light_data, 2));
            for c = 1:size(obj.light_data, 1)
                light = rescale(obj.light_data(c, :));
                
                
                dark = rescale(obj.dark_data(c, :));
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
                pause
            end
        end
        
        function out = getClusters(obj)
            out = obj.clust_id;
        end
        
        function clustID = cluster(obj, data)
            % From Michael 20Dec2019
            
            varExpl = 0.98;    % criteria for variance explained (deteremines # of PCs)
            
            numSamp = size(data, 2);
            numCells = size(data, 1);
            
            %% PCA
            [~,score,latent] = pca(data);
            PC_var = cumsum(latent/sum(latent));
            numPCs = find(PC_var > varExpl,1);
            
            %% Cluster responses
            PCA_resp = score(:,1:numPCs);
            Z = linkage(PCA_resp,'ward','euclidean');
            dendrogram(Z)
            set(gcf,'color',[1 1 1])
            numClust = input('How many clusters? ');    % Determine cluster number on dendrogram
            %saveas(gcf,'Cluster_dendrogram','fig')
            close all
            clustID = cluster(Z,'maxclust',numClust);
            
            %% Recluster by correlation to mean response (to clean up sorting)
            clusterMeanTraces = zeros(numClust, numSamp);
            for c = 1:numClust % calculate average traces
                clusterMeanTraces(c,:) = mean(data(clustID==c,:));
            end
            for n = 1:numCells % assign ased on CC to average traces
                curr_trace = data(n,:);
                for c = 1:numClust
                    CC = corrcoef(curr_trace,clusterMeanTraces(c,:));
                    clustCC(c) = CC(2);
                end
                [~,clustID(n)] = max(clustCC);
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
            
            figure(1)
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
    end
end