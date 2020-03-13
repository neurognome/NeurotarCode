classdef HeadDirectionAnalyzer < handle

    % Make sure you have the circular statistics toolbox
    properties (Constant = true)
        RADIUS    double = 125; % Radius of the arena
    end
    
    properties
        data struct
        floating struct
        mean_quadrant_correlation double
        direction_info double
        bilobality_fit_info struct
        %   binned_data double
        
        is_head_direction logical
        is_bilobed logical
        is_unilobed logical
        
        fix_heading_flag logical = false;
        bin_width double = 5; % previously 3
        shuffled_threshold double
    end
    
    properties (Access = public)
        init_data struct
        init_floating struct
        is_moving logical
        remove_slow_flag logical = false;
    end
    
    methods % Main methods
        function obj = HeadDirectionAnalyzer(data, floating)
            obj.init_data = data;
            obj.init_floating = floating;
            
            obj.data = data;
            obj.floating = floating;
            
            obj.extractHeading();
        end

        function printStats(obj)
            fprintf('%0.2f%% cells tuned\n', mean(obj.is_head_direction) * 100);
        end
        
        function out = calculatePreferredDirection(obj, method, binned_data)
            % Get the preferred direction of each cell, either vectorsum or just the max
            if nargin < 2 || isempty(method)
                method = 'vectorsum';
                fprintf('No method provided, defaulting to vectorsum\n')
            end
            
            if nargin < 3 || isempty(binned_data)
                binned_data = obj.binData([], [], [], true);
            end
            
            out = zeros(size(binned_data, 1), 2);
            switch method
            case 'vectorsum'
                for c = 1:size(binned_data, 1)
                    dat = binned_data(c, :);
                        dat(dat < 0) = 0; % rectify
                     %   dat = rescale(dat);
                        %dat(dat < 0) = 0; % rectify
                        out(c, :) = obj.calculateVectorSum(dat, true); % True here because we want the wrapped one
                    end
                    out(:, 1) = wrapTo2Pi(out(:, 1));
                case 'max'
                    for c = 1:size(binned_data, 1)
                        [out(c, 2) ,out(c, 1)] = max(binned_data(c, :));
                    end
                otherwise
                    fprintf('Invalid method provided (methods: vectorsum or max)\n')
                    return
                end

                obj.direction_info = out;
            end

            function calculateHeadDirectionIdx_direction(obj)
            % The preferred direction here is calculated by splitting the recording into sections (generally quadrants), then
            % seeing how "stable" the responses are across the sections. A correlation greater than 0.2 is considered to be
            % "tuned"
            % the following analyses are based on Giocomo et al 2017, in
            % current biology...
            
            neural_data = obj.initializeData(obj.data,'spikes');
            
            obj.direction_info = obj.calculatePreferredDirection('vectorsum');
            
            obj.is_head_direction = obj.direction_info(:, 2)' > 0.2;% temporary thershold obj.shuffled_threshold;
        end

        function calculateHeadDirection(obj, iterations, threshold)
            if nargin < 2 || isempty(iterations)
                iterations = 100;
            end

            if nargin < 3 || isempty(threshold)
                threshold = 95;
            end

            data = obj.binData();
            [~, ori_pref] = max(data, [], 2);
            osi_vec = obj.calculateOSI(data, ori_pref);

            % shuffle 
            for i = 1:iterations
                shuffled_data = circshift(data, randi(size(data, 2)), 2); % rotate
                shuffled_osi(:, i) = obj.calculateOSI(shuffled_data, ori_pref);
            end
            threshold = prctile(shuffled_osi, 95, 2);
            obj.is_head_direction = osi_vec > threshold';
        end

        function osiVec = calculateOSI(obj, binned_data, ori_pref)                
                % Prepare data
               % oriPref = obj.getPreferredDirection();
              %  oriPref = round(oriPref .* size(binned_data, 2) ./ (2*pi));
              n_ori = size(binned_data, 2);
              if nargin < 3 || isempty(ori_pref)
                [~, ori_pref] = max(binned_data,[],2);
            end

            ct = 0;
            win = 2;
            for resp = binned_data'
                ct = ct + 1;
                pref = mean(resp(max([1, ori_pref(ct) - win]) : min([length(resp), ori_pref(ct) + win])));
                orth = mean(resp([max([1, mod(ori_pref(ct) - n_ori / 4 - 1, n_ori) - win + 1]):...
                    min([length(resp), mod(ori_pref(ct) - n_ori / 4 - 1, n_ori) + win + 1]),...
                    max([1, mod(ori_pref(ct) + n_ori / 4 - 1, n_ori) - win + 1]):...
                    min([length(resp), mod(ori_pref(ct) + n_ori / 4 - 1, n_ori) + win + 1])]));

                osiVec(ct) = (max([1 pref])-max([1 orth]))/(max([1 pref])+max([1 orth]));
            end
        end
        
        function removeMovingSamples(obj, speed_threshold, peak_width, search_window)
            % This code is used to remove moving samples from the data,  here just defined as speed less than 10, it's pretty
            % lenient
            
            if nargin < 2
                speed_threshold = 5; % Default
            end
            if nargin < 3
                peak_width = 5;
            end
            if nargin < 4
                search_window = 20;
            end
            
            obj.findMovingTimes(speed_threshold, peak_width, search_window);
            obj.remove_slow_flag = true;
        end
    end
    
    methods % Beta analyses
        function flip_score = calculateFlipScore(obj)
            % Get the binned data and calculate autocorrelation
            binned_data = obj.getBinnedData();
            % If we rotate by 6 degree increments, that means to move 2 bins each time
            rotations = 0:2:size(binned_data, 2);
            rotations = rotations(1:end-1); % Don't need to overlap the last one twice ya feel?
            
            auto_corr = zeros(size(binned_data, 1), length(rotations));
            for i_cell = 1:size(binned_data, 1)
                data = binned_data(i_cell, :);
                for i_rot = 1:length(rotations)
                    auto_corr(i_cell, i_rot) = corr(circshift(data, rotations(i_rot))', data');
                end 
            end
            
            % Maybe we should add a little "leeway" here
            flip_score = auto_corr(:, floor((180/360) * length(rotations)) + 1) - ...
            auto_corr(:, floor((90/360) * length(rotations)) + 1);
            
            %The following analyses were using "standard" autocorrelational measures, but we're going to try the method
            %outlined in the paper
            %{ 
            for i_cell = 1:size(binned_data, 1)
                % auto_corr(i_cell, :) = xcorr(binned_data(i_cell, :), binned_data(i_cell, :), 'coeff');
                auto_corr(i_cell, :) = ifft(fft(binned_data(i_cell, :)).*conj(fft(binned_data(i_cell, :)))); % Calculate the
                % autocorrelation
            end
            
            n_bins = size(binned_data, 2);
            auto_corr = auto_corr ./ auto_corr(:, n_bins/2 + 1); % Normalized so that 0 lag = 1;
            
            % Calculate flip score
            flip_score = auto_corr(:, (180 / 360) * n_bins + 1) - auto_corr(:, (90 / 360) * n_bins + 1);
            
            flip_score(abs(flip_score) > 5) = 0;  % Almost always due to errors
            %}
        end
    end
    

    methods % Setters and getters
        function out = getHeading(obj)
            out = obj.floating.heading;
        end
        
        function out = getBinnedData(obj)
            out = obj.binData();  % Call w/o input arguments so it correctly gets the thing
        end
        
        function obj = setHeadingFlag(obj, val)
            obj.fix_heading_flag = val;
            obj.getHeading();
        end
        
        function obj = setBinWidth(obj, bin_width)
            obj.bin_width = bin_width;
        end
        
        function out = getPreferredDirection(obj)
            out = obj.direction_info(:, 1);
        end 
        
        function out = getTimeSeries(obj)
            out = obj.data.spikes;
        end
    end
    
    methods (Access = protected) % Helper methods

        function [x, y] = getCoords(obj)
            % Convetrs the raw X and Y from floating into usable coordinates for plotting
            x = round(obj.floating.X + obj.RADIUS)';  % Round to turn into coordinates
            y = round(obj.floating.Y + obj.RADIUS)';
            
            % Bound at radius
            x = max(x, 1);
            x = min(x, 2 * obj.RADIUS);
            
            y = max(y, 1);
            y = min(y, 2 * obj.RADIUS);
        end
        
        function extractHeading(obj)
            % Converts alpha to heading via some transformations we figured out
            if obj.fix_heading_flag
                obj.floating.heading = calculateHeading(...
                    obj.floating.r, obj.floating.phi, obj.floating.alpha, obj.RADIUS);
            else
                obj.floating.heading = obj.floating.alpha;
            end
            
            function h = calculateHeading(r, phi, alpha, radius)
                h = 2 .* asind((r .* sind(alpha - phi)) ./ (2 .* radius)) + alpha;
                h = wrapTo180(h); % to account for passing 180 on the thing
            end
        end
        
        function varargout = initializeData(obj, data_struct, varargin)
            % Initializes the data in preparation for multiple analyses. Mainly checks if you want to remove moving samples
            % or not and goes from there
            for ii = 1:length(varargin)
                d = data_struct.(varargin{ii});
                
                if size(d,1) > size(d, 2)
                    d = d';
                end
                
                if obj.remove_slow_flag
                    varargout{ii} = d(:, obj.is_moving);
                else
                    varargout{ii} = d;
                end
            end
        end
        
        function findMovingTimes(obj, speed_threshold, peak_width, search_window)
            % Finds times that are considered moving based on a simple thershold, and reports them
            if ~isempty(obj.is_moving)
                disp('Already removed moving times, don''t run it again')
                return
            end
            
            speed = obj.initializeData(obj.floating, 'speed');
            obj.is_moving = speed > speed_threshold;
            moving_idx = find(obj.is_moving);
            init_moving = obj.is_moving; % store init
            for idx = moving_idx
                pre_window = init_moving(max(1, (idx - search_window)):idx-1);
                post_window = init_moving((idx+1):min(length(init_moving), idx + search_window));
                obj.is_moving(idx) = ~((sum(pre_window) < peak_width && sum(post_window) < peak_width));
            end
        end
        
        function out = binData(obj, bin_width, data, heading, fold_flag)
            % Bins data based on a determined bindwidth
            
            % added a bunch of stuff in case we want folding.. this is b/c we only want to fold sometimes
            if nargin < 2 || isempty(bin_width)
                bin_width = obj.bin_width;
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
                    temp = data(c, groups == u_groups(bin));
                    % histogram(temp)
                    % pause
                    out(c, bin) = mean(temp); % should I only take values that are nonzero?
                    % out(c, bin) = mean(data(c, groups  == u_groups(bin)));
                end
            end
            
            out(isnan(out)) = 0; % NanCheck
            
            % wrap to 2pi
            if fold_flag
                out = mean(cat(3, out(:, 1:length(u_groups)/2), out(:, length(u_groups)/2 + 1:end)), 3);
            end
            out = movmean(out, 15/obj.bin_width, 2); % 10 bins, 5 on each side, = 15 degree on each side, same as Giocomo et al 2014
            
        end
    end
    
    methods (Static = true) % Static methods
        function out = calculateVectorSum(data, fold_flag)
            if nargin < 2 || isempty(fold_flag)
                fold_flag = false; % By default, let's fold it to account for bipolar neurons as well
            end
            
            fold_flag = false; % moved the angle doubling elsewhere, but keeping here for now just in case

            getHorz = @(v, theta) v .* cos(theta);
            getVert = @(v, theta) v .* sin(theta);
            getAng = @(vert, horz) atan2(vert, horz);
            getMag = @(vert, horz) sqrt(horz ^ 2 + vert ^ 2);
            
            if fold_flag
                data = mean([data(1:size(data, 2)/2); data(size(data, 2)/2 + 1:end)]);
                theta = linspace(0, 2*pi, length(data));
                warning('Because of the angle doubling, don''t particularly trust the actual values for pref dir')
            else
                theta = linspace(0, 2*pi, length(data));
            end
            
            %  data = rescale(data); % to rescale or not to rescale....
            
            
            % testing data : dat = circshift([9 8 7 6 5 6 7 8],ii);
            h = getHorz(data, theta);
            v = getVert(data, theta);
            
            % Changed from sum to mean, shouldn't change anything... more similar to Giocomo
            r_h = mean(h);
            r_v = mean(v);
            
            m = getMag(r_v, r_h);
            a = getAng(r_v, r_h);
            
            %% Correcting for quadrant errors
            if r_h > 0
                ang = a;
            elseif r_h < 0
                ang = a + pi;
            else
                disp('Uh... that''s not supposed to happen')
                ang = 0;
            end
            
            out = [ang, m];
        end
        
        function distance = getCircularDistance(loc1, loc2)
            % ensure p2 >= p1 always
            if (loc2 < loc1)
                t = loc1(1); loc1 = loc2;  loc2 = t;
            end
            % forward dist
            dist = loc2 - loc1 + 1;
            % rev dist
            dist2 = loc1 - 1 + 12 - loc2;
            distance = min(dist, dist2);
        end
    end
end


%{
Deprecated code:
        function quadrantAnalysisHeadDirection(obj, plot_flag)
            % The purpose of this method is to use the preferred direction of the cell to figure out what other preferred
            % directions it may have. Essentially, for a bilobed tuning cell, we only expected bilobality in certain
            % directions
            
            if nargin < 2
                plot_flag = false;
            end
            
            binned_data = obj.binData();
            obj.calculatePreferredDirection('vectorsum');  % Just to make sure we actually have vectorsum
            
            north_L = [3*pi/2 + pi/4, 2*pi]; % North needs to split because of cyclic nature
            north_R = [0, pi/2 - pi/4];
            east = [pi/2 - pi/4, pi/2 + pi/4];
            south = [pi - pi/4, pi + pi/4];
            west = [3*pi/2 - pi/4, 3*pi/2 + pi/4];
            
            pref_dir = obj.direction_info(:, 1);
            
            north_pref =  (pref_dir > north_L(1) & pref_dir <= north_L(2)) | ...
                (pref_dir >= north_R(1) & pref_dir < north_R(2));
            east_pref = pref_dir >= east(1) & pref_dir < east(2);
            south_pref = pref_dir >= south(1) & pref_dir < south(2);
            west_pref = pref_dir >= west(1) & pref_dir < west(2);
            
            if plot_flag
                % Visualize
                figure
                
                subplot(2, 4, 1)
                imagesc(binned_data(north_pref & obj.is_head_direction, :))
                title('North')
                axis image
                subplot(2, 4, 5)
                plot(mean(binned_data(north_pref & obj.is_head_direction, :)));
                prettyPlot()
                
                subplot(2, 4, 2)
                imagesc(binned_data(east_pref & obj.is_head_direction, :))
                title('East')
                axis image
                subplot(2, 4, 6)
                plot(mean(binned_data(east_pref & obj.is_head_direction, :)));
                prettyPlot()
                
                subplot(2, 4, 3)
                imagesc(binned_data(south_pref & obj.is_head_direction, :))
                title('South')
                axis image
                subplot(2, 4, 7)
                plot(mean(binned_data(south_pref & obj.is_head_direction, :)));
                prettyPlot()
                
                subplot(2, 4, 4)
                imagesc(binned_data(west_pref & obj.is_head_direction, :))
                title('West')
                axis image
                subplot(2, 4, 8)
                plot(mean(binned_data(west_pref & obj.is_head_direction, :)));
                prettyPlot()
            end
        end

        
        function trajectoryActivityPlot(obj, c)
            % Used to check the activity of the given cell at each position. Placecell-type stuff
            [x, y] = obj.getCoords();
            z = zeros(size(x));
            col = obj.data.spikes(c, :);
            col(col < 1) = min(col(:));
            surface([x;x],[y;y],[z;z],[col;col],...
                'AlphaData', [col;col],...
                'EdgeAlpha', 'interp',...
                'facecol','no',...
                'edgecol','interp',...
                'linew',2);
            colormap jet
            axis square
            title(sprintf('Cell %d', c))
        end
        
        function pointActivityPlot(obj, c)
            % Essentially a thresholded version of trajectoryActivityPlot.
            [x, y] = obj.getCoords();
            high_activity_frames = find(obj.data.spikes(c, :) > prctile(obj.data.spikes(c, :), 99));
            
            arena = zeros(2 * obj.RADIUS, 2 * obj.RADIUS);
            for f = high_activity_frames
                arena(x(f), y(f)) = 1;
            end
            
            imagesc(imgaussfilt(arena, 0.5))
            set(gca, 'YDir', 'normal')
            axis square
            colormap hot
        end

        
        function calculateHeadDirectionIdx(obj, n_sections)
            % The preferred direction here is calculated by splitting the recording into sections (generally quadrants), then
            % seeing how "stable" the responses are across the sections. A correlation greater than 0.2 is considered to be
            % "tuned"
            % the following analyses are based on Giacomo et al 2017, in
            % current biology...
            if nargin < 2
                n_sections = 4;
            end
            
            neural_data = obj.initializeData(obj.data,'spikes');
            heading = obj.initializeData(obj.floating,'heading');
            
            quadrant_correlations = obj.calculateQuadrantCorrelations(neural_data, heading, n_sections);
            
            obj.mean_quadrant_correlation = mean(quadrant_correlations,2);
            
            if isempty(obj.shuffled_threshold)
                obj.getShuffledThreshold(neural_data, heading, n_sections);
            end
            
            obj.is_head_direction = obj.mean_quadrant_correlation' > obj.shuffled_threshold;
        end

  
        function [coeffs, r_sqr] = temporaryMegaFit(obj, data)
            bound = @(x, bl, bu) min(max(x, bl), bu);
            %
            modelfun = @(b, x) b(1) + ...
                max([0, b(2)]) * exp(-(x(:, 1) - bound(b(3), 0, x(end, 1))) .^ 2 / (2 * b(4) .^ 2)) + ...
                max([0, b(5)]) * exp(-(x(:, 1) - bound(b(6), 0, x(end, 1))) .^ 2 / (2 * b(4) .^ 2)); % Based on Michael natneuro paper
            
            % fit the function
            for ii = 1:size(data, 3)
                textprogressbar(ii / size(data, 3))
                [coeffs{ii}, r_sqr{ii}] = obj.fitDoubleGaussian(data(:, :, ii), modelfun);
            end
        end
        
        function out = bilobalityAnalysis(obj, check_flag)
            % Here we're going to find some kind of metric to measure bilobality.
            
            if nargin < 2
                check_flag = false;
            end
            
            binned_data = obj.binData();
            bound = @(x, bl, bu) min(max(x, bl), bu);
            %             modelfun = @(b, x) b(1) + ...
            %                 max([0, b(2)]) * exp(-(x(:, 1) - bound(b(3), 0, x(end, 1))) .^ 2 / (2 * b(4) .^ 2)) + ...
            %                 max([0, b(5)]) * exp(-(x(:, 1) - bound(b(6), 0, x(end, 1))) .^ 2 / (2 * b(4) .^ 2)); % Based on Michael natneuro paper
            modelfun = @(b, x) b(1) + ...
                max([0, b(2)]) * exp(-(x(:, 1) - b(3)) .^ 2 / (2 * b(4) .^ 2)) + ...
                max([0, b(5)]) * exp(-(x(:, 1) - b(6)) .^ 2 / (2 * b(4) .^ 2)); % Based on Michael natneuro paper
            
            % fit the function
            [coeffs, r_sqr] = obj.fitDoubleGaussian(binned_data, modelfun);
            
            % What criteria do we want?
            %   1) There needs to be a minimum distance between the gaussians
            %       Accounted for in the model
            %   2) Well fit by a 2 gaussian model
            %       Accounted for by r squared metric
            %   3) Existence of two lobes, single lobed responses stored elsewhere
            
            % Adjusting to get rid of values lower than the base possible (fitting errors)
            %             coeffs(:, 3) = max(coeffs(:, 3), 1);
            %             coeffs(:, 6) = max(coeffs(:, 6), 1);
            
            % Going to have to play with these parameters to get it to work right...
            for i_coeff = 1:size(coeffs, 1)
                peak_distances(i_coeff) = obj.getCircularDistance(coeffs(i_coeff, 3), coeffs(i_coeff, 6));
            end
            
            
            % think about circular distance
            is_separated = (peak_distances > 3)';
            is_well_fit = r_sqr > 0.5;
            is_double_lobed = coeffs(:, 5) > 0.1 * coeffs(:, 2); % Second peak needs to be at least 20% of the first peak
            obj.is_bilobed = is_well_fit & is_double_lobed & is_separated;
            obj.is_unilobed = is_well_fit & ~is_double_lobed;
            
            
            % obj.is_bilobed = is_well_fit;
            out.coeffs = coeffs;
            out.r_squared = r_sqr;
            obj.bilobality_fit_info = out;
            
            bilob_idx = coeffs(obj.is_bilobed & obj.is_head_direction, 3);
            %[~, bilob_idx] = max(bilob_data, [], 2);
            tabulate(round(bilob_idx))
            
            
            if check_flag
                x_vals = linspace(1, size(binned_data, 2), 100);
                for ii = 1:size(binned_data, 1)
                    if (obj.is_bilobed(ii) || obj.is_unilobed(ii)) && obj.is_head_direction(ii)
                        y = modelfun(coeffs(ii, :), x_vals');
                        plot(x_vals, y)
                        hold on
                        plot(binned_data(ii, :))
                        title([num2str(ii), ' | ' num2str(r_sqr(ii)), 'Bilobed: ', num2str(obj.is_bilobed(ii)), ...
                            'Unilobed: ', num2str(obj.is_unilobed(ii))])
                        pause
                        hold off
                    end
                end
            end
        end


function saveAnalysisData(obj)
            mean_quadrant_correlation = obj.mean_quadrant_correlation;
            direction_info = obj.direction_info;
            is_head_direction = obj.is_head_direction;
            save_filename = sprintf('HDA_data_%s.mat', date);
            save(save_filename, 'mean_quadrant_correlation',...
                'direction_info',...
                'is_head_direction')
            fprintf('File saved as: %s\n', save_filename);
        end
         methods % Visualization methods
        function bilobalityPlot(obj)
            binned_data = obj.getBinnedData();
            bilob_idx = obj.bilobality_fit_info.coeffs(obj.is_bilobed & obj.is_head_direction, 3);
            bilob_data = binned_data(obj.is_bilobed & obj.is_head_direction, :);
            [r_max, max_idx] = max(bilob_data, [], 2);
            r_min = min(bilob_data, [], 2);
            
            plot_data = rescale(bilob_data, 'InputMin', r_min, 'InputMax', r_max);
            
            for d = 1:size(plot_data, 1)
                plot_data(d, :) = circshift(plot_data(d, :), 1 - max_idx(d));
            end
            
            figure
            subplot(2, 1, 1)
            imagesc(plot_data);
            ylabel('Cell #')
            xlabel('Direction')
            
            subplot(2, 1, 2)
            %theta = (bilob_idx / length(bilob_idx)) .* 2 * pi;
            histogram(bilob_idx, [0.5:12.5]);  % Force the bin edges so that it lines up properly
            axis([0 13 ylim])
            ylabel('Counts')
            xlabel('Bilobality Direction')
            tabulate(round(bilob_idx))
        end
        
        function alignedPlot(obj)
            obj.getPreferredDirection('max');
            % Aligns to preferred direction, and then plots them, to check for second peaks mainly
            pref_dir = obj.direction_info(:, 1); % prefdir
            binned_data = binned_data(obj.mean_quadrant_correlation > 0.2, :);
            pref_dir = pref_dir(obj.mean_quadrant_correlation > 0.2);
            aligned_data = zeros(size(binned_data));
            for c = 1:size(binned_data, 1)
                aligned_data(c, :) = rescale(circshift(binned_data(c, :), 1 - pref_dir(c)));
                
                temp = findpeaks(aligned_data(c, :), 'SortStr', 'descend');
                bimodality_idx(c) = aligned_data(c, 1) / temp(1);
            end
            [~, sortingIdx] = sort(bimodality_idx);
            
            aligned_data = aligned_data(sortingIdx, :);
            
            % sort
            figure;
            subplot(1, 2, 1)
            imagesc(aligned_data)
            subplot(1, 2, 2)
            plot(mean(aligned_data));
        end
        
        function polarPlot(obj, c)
            % Simple polar-plotter for single cells
            plt = RawDataPlots();
            plt.setData(obj.binData());
            plt.polarPlot(c, 'LineWidth', 5)
            title(sprintf('Cell %d', c))
        end
    end
            function out = getSegmentedData(obj, n_sections)
            neural_data = obj.initializeData(obj.data,'spikes');
            heading = obj.initializeData(obj.floating,'heading');
            segmented_data = obj.segmentData(neural_data, n_sections);
            segmented_heading = obj.segmentData(heading, n_sections);
            
            out = obj.binSegmentedData(segmented_data, segmented_heading, n_sections);    
        end

        function out = getSignificantCells(obj, threshold)
            if nargin < 2
                threshold = 0.2;  % From the paper
            end
            out = find(obj.mean_quadrant_correlation > threshold)';
        end
                function q_binned_data = binSegmentedData(obj, segmented_data, segmented_heading, n_sections)
            q_binned_data = zeros(size(segmented_data{1}, 1), length(-180:obj.bin_width:180) - 1, n_sections);
            for q = 1:n_sections
                q_binned_data(:, :, q) = obj.binData(obj.bin_width, segmented_data{q}, segmented_heading{q});
            end
        end
        

        function getShuffledThreshold(obj, data, threshold)
            if nargin < 3 || isempty(threshold)
                threshold = 99;
            end
            
            disp('Shuffling data...')
            
            for iter = 1:100 % 100 iterations
                shuffled_neural_data = obj.shuffleData(data);
                shuffled_binned_data = obj.binData([], shuffled_neural_data);
               % shuffled_quadrant_correlations(:, :, iter) = ...
                %    obj.calculateQuadrantCorrelations(shuffled_neural_data, heading, n_sections);
                shuffled_direction_info(:, :, iter) = obj.calculatePreferredDirection('vectorsum', shuffled_binned_data);       
            end
            obj.shuffled_threshold = prctile(squeeze(shuffled_direction_info(:, 2, :))', threshold);
        end
   
        
        function quadrant_correlations = calculateQuadrantCorrelations(obj, neural_data, heading, n_sections)
            q_data = obj.segmentData(neural_data, n_sections);
            q_heading = obj.segmentData(heading, n_sections);     
            q_binned_data = obj.binSegmentedData(q_data, q_heading, n_sections);
            % compare quadrants in pairwise (1-1,1-2,1-3,1-4,2-3,2-4,3-4)
            possible_combinations = nchoosek(1:n_sections, 2);
            quadrant_correlations = zeros(size(neural_data, 1), size(possible_combinations, 1));
            for ii = 1:size(possible_combinations, 1)
                for c = 1:size(neural_data, 1)
                    quadrant_correlations(c, ii) = corr(...
                        q_binned_data(c, :, possible_combinations(ii, 1))',...
                        q_binned_data(c, :, possible_combinations(ii, 2))');
                end
            end
        end     

                
        function segmented_data = segmentData(data, n_sections)
            total_dur = size(data, 2);
            section_length = round(total_dur./n_sections);
            
            segmented_data = cell(n_sections, 1);
            for q = 1:n_sections % Separate timeseries into sections
                segmented_data{q}    = data(:, ...
                    (q - 1) * section_length + 1 : min(q*section_length, size(data, 2)));
            end
        end

function shuffled_data = shuffleData(data)
            % This function shuffles the given data. Please provide data in cells x data pts
            shuffled_data = zeros(size(data));
            for row = 1:size(data, 1)
                shift_amt = randi([20, size(data, 2) - 20]);
                shuffled_data(row, :) = circshift(data(row, :), shift_amt);
            end
        end
%}