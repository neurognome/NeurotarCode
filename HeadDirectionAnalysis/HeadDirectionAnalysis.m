classdef HeadDirectionAnalysis < RawDataPlots & SummaryDataPlots
    
    % Make sure you have the circular statistics toolbox
    properties (Constant = true)
        radius    double = 120; % Radius of the arena
    end
    
    properties
        DFF double
        heading double
        numCells double
        
        spike_flag logical
        neural_data
        
        fixHeadingFlag logical = 1;
        thresh    double = 0.05; % p-value threshold for significant tuning
        binWidth double = 20;
    end
    
    properties (Access = protected)
        initData struct
        initFloating struct
    end
    
    methods
        
        function obj = HeadDirectionAnalysis(data, floating)
            obj.initData = data;
            obj.initFloating = floating;
            [obj.neural_data, obj.heading] = obj.initializeData(data, floating);
        end
        
        function findHeadDirectionCells(obj) % just a wrapper for full analysis later
        end
       
        function calculateHeadDirectionIdx(obj, numSections)
            % the following analyses are based on Giacomo et al 2017, in
            % current biology...
            if nargin < 2
                numSections = 4;
            end
            
            totalDur = length(obj.heading);
            quadrantLength = round(totalDur./numSections);
            
            q_data = cell(numSections, 1);
            q_heading = cell(numSections, 1);
            for q = 1:4 % 4 quadrants...
                if q == 4 % because it might not equally divide into 4, we might need to give one a lil' extra
                    q_data{q} = obj.data(:, (q - 1) * quadrantLength + 1 : end);
                    q_heading{q} = obj.heading((q - 1) * quadrantLength + 1 : end);
                else
                    q_data{q}    = obj.data(:,(q - 1) * quadrantLength + 1 : q*quadrantLength);
                    q_heading{q} = obj.heading((q - 1) * quadrantLength + 1 : q*quadrantLength);
                end
            end
            
            
            q_binnedDFF = zeros(size(obj.binDFF, 1), size(obj.binDFF, 2), numSections);
            % Calculate the head direction preference per quadrant
            for q = 1:numSections
                q_binnedDFF(:, :, q) = obj.binDFF(q_data{q}, q_heading{q});
            end
            
            % compare quadrants in pairwise (1-1,1-2,1-3,1-4,2-3,2-4,3-4)
            possibleCombinations = nchoosek(1:4, 2);
            
            quadrantCorrelations = zeros(size(obj.data, 1), size(possibleCombinations, 1));
            for c = 1:size(obj.data, 1)
                for ii = 1:size(possibleCombinations, 1)
                    quadrantCorrelations(c, ii) = corr(q_binnedDFF(c, :, possibleCombinations(ii,1))', q_binnedDFF(c, :, possibleCombinations(ii, 2))');
                end
            end

        end
        
        function shuffleHeadDirectionIdx(obj, data, floating)
            
            % where should we shuffle these data?
        end
        %% Setters and Getters
        function setHeadingFlag(obj, val)
            obj.fixHeadingFlag = val;
            [obj.DFF, obj.heading] = obj.initializeData(obj.initData, obj.initFloating);
        end
    end
    
    methods (Access = public)

        function [neural_data,heading] = initializeData(obj,data,floating)
            try % temporarily to account of spikeInferenc/noSpikeInference
                neural_data = data.spikes;
                obj.spike_flag = true;
            catch
                fprintf('Not spike inferred, using DFF\n')
                neural_data = data.DFF;
                obj.spike_flag = false;
            end
            X = floating.X;
            Y = floating.Y;
            alpha = floating.alpha;
            phi = floating.phi;
            r = floating.r;
            
            if obj.fixHeadingFlag
                heading = obj.calculateHeading(r,phi,alpha);
            else
                heading = alpha;
            end
        end
        
        function [data, floating] = removeStill(obj, speed_threshold)
            if nargin < 2
                speed_threshold = 10; % Default
            end
            
            isTooSlow = obj.initFloating.speed < speed_threshold;

            floating_fields = fields(floating); % Get fields of floating

            for ii = 1:length(floating_fields) % Cut it 
                if ~strcmp(floating_fields{ii},'time') % don't mess with time
                floating.(floating_fields{ii})(isTooSlow) = [];
                end
            end
            
            if obj.spike_flag
                data.spikes(:,isTooSlow) = [];
            else
                data.DFF(:,isTooSlow) = [];
            end
            
            [obj.neural_data,obj.heading] = obj.initializeData(data,floating);
        end
        
        function out = detectCells(obj)
            for ii = 1:size(obj.data, 1)
                otest(ii) = circ_otest(obj.data(ii, :)); %nonuniformity
                if otest(ii) < obj.thresh
                    rtest(ii) = circ_rtest(obj.data(ii, :)); %unimodality
                else
                    rtest(ii) = 1;
                end
            end
            out = (otest<obj.thresh) & (rtest<obj.thresh);
        end
        
        function out = getPreferredDirection(obj) % We should consider using the vector sum here too
            fprintf('Getting preferred directions \n');
            for ii = 1:size(obj.data, 1)
                % lm = fit([1:length(data(ii,:))]',data(ii,:)','gauss1','Upper',[Inf 36 Inf], 'Lower', [-Inf 0 -Inf]);
                %  out(ii) = lm.b1;
                [~,out(ii)] = max(obj.data(ii, :));
            end
            out = out * obj.binWidth - obj.binWidth; % convert to degrees, accounting for the extra 10 you get because -180:10:180
        end
        
        function out = binDFF(obj, data, heading)
            
            if nargin < 3
                data = obj.data;
                heading = obj.heading;
            end
            binEdges = -180:obj.binWidth:180; % Because alpha is [-180,180];
            groups = discretize(heading, binEdges);
            u_groups = unique(groups);
            for c = 1:size(data, 1)
                for bin = 1:length(u_groups)
                    out(c, bin) = mean(data(c, groups == u_groups(bin)));
                end
            end
        end
        
        function h = calculateHeading(obj,r,phi,alpha)
            h = 2.*asind((r.*sind(alpha - phi))./(2.*obj.radius)) + alpha;
            h = wrapTo180(h); % to account for passing 180 on the thing
        end
        
    end
    
end



%% Major analysis code

% when you're done messing around in the script, you can incorporate as a big thing, but until then, run the individual functions in the script ok

%{
function findHeadDirectionCells(obj)
            
            binned_DFF = obj.binDFF(obj.workingData.get('heading'),obj.workingData.get('DFF'));
            obj.analysisData.add(binned_DFF);
            
            p = obj.detectCells(obj.workingData.get('DFF'));
            pref_dir = obj.getPreferredDirection(obj.analysisData.get('binned_DFF'));
            
            obj.analysisData.add('p','pref_dir');
            
            isDirectionTuned = p < obj.thresh;
            obj.analysisData.add(isDirectionTuned);

        end
%}

