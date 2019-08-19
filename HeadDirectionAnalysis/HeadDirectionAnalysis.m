classdef HeadDirectionAnalysis < handle
    
    % Make sure you have the circular statistics toolbox
    properties (Constant = true)
        radius    double = 120; % Radius of the arena
    end
    
    properties
        data struct
        floating struct
        meanQuadrantCorrelation double
       
        fixHeadingFlag logical = false;
        binWidth double = 20;
    end
    
    properties (Access = public)
        initData struct
        initFloating struct
        isMoving logical
        removeSlowFlag logical = false;
    end
    
    methods
        function obj = HeadDirectionAnalysis(data, floating)
            obj.initData     = data;
            obj.initFloating = floating;
            
            obj.data     = data;
            obj.floating = floating;
            
            obj.get_heading();
        end
        
        function calculate_head_direction_idx(obj, numSections)
            % the following analyses are based on Giacomo et al 2017, in
            % current biology...
            if nargin < 2
                numSections = 4;
            end
            
            neural_data = obj.initializeData(obj.data,'spikes');
            heading     = obj.initializeData(obj.floating,'heading');
            
            totalDur = length(heading);
            sectionLength = round(totalDur./numSections);
            
            q_data = cell(numSections, 1);
            q_heading = cell(numSections, 1);
            for q = 1:numSections % 4 quadrants...
                    q_data{q}    = neural_data(...
                            :,(q - 1) * sectionLength + 1 : min(q*sectionLength, length(neural_data)));
                    q_heading{q} = heading(...
                            (q - 1) * sectionLength + 1 : min(q*sectionLength, length(heading)));
            end
            q_binnedDFF = zeros(size(neural_data, 1), length(-180:obj.binWidth:180)-1, numSections);
            % Calculate the head direction preference per quadrant
            for q = 1:numSections
                q_binnedDFF(:, :, q) = obj.bin_DFF(q_data{q}, q_heading{q});
            end
            % compare quadrants in pairwise (1-1,1-2,1-3,1-4,2-3,2-4,3-4)
            possibleCombinations = nchoosek(1:4, 2);
            quadrantCorrelations = zeros(size(neural_data, 1), size(possibleCombinations, 1));
            for ii = 1:size(possibleCombinations, 1)
                for c = 1:size(neural_data, 1)
                    quadrantCorrelations(c, ii) = corr(...
                        q_binnedDFF(c, :, possibleCombinations(ii,1))',...
                        q_binnedDFF(c, :, possibleCombinations(ii, 2))');
                end
            end
            obj.meanQuadrantCorrelation = mean(quadrantCorrelations,2);
        end
        
        function find_moving_samples(obj, speed_threshold, peakWidth, searchWindow)
            if nargin < 2
                speed_threshold = 10; % Default
            end
            if nargin < 3
                peakWidth = 5;
            end
            if nargin < 4
                searchWindow = 10;
            end
            
            obj.find_moving_times(speed_threshold, peakWidth, searchWindow);
            obj.removeSlowFlag = true;
        end
             
        %% Setters and Getters
        function set_heading_flag(obj, val)
            obj.fixHeadingFlag = val;
            obj.get_heading();
        end
    end
    
    methods (Access = protected)
        function get_heading(obj)
            if obj.fixHeadingFlag
                obj.floating.heading = obj.calculate_heading(...
                    obj.floating.r, obj.floating.phi, obj.floating.alpha);
            else
                obj.floating.heading = obj.floating.alpha;
            end
        end
        
        function varargout = initializeData(obj, dataStruct, varargin)
            for ii = 1:length(varargin)
                d = dataStruct.(varargin{ii});
                
                if size(d,1) > size(d, 2)
                    d = d';
                end
                
                if obj.removeSlowFlag
                    varargout{ii} = d(:, obj.isMoving);
                else
                    varargout{ii} = d;
                end
            end
        end
        
        function find_moving_times(obj, speed_threshold, peakWidth, searchWindow)
            speed = obj.initializeData(obj.floating, 'speed');
            obj.isMoving = speed > speed_threshold;    
            movingIdx = find(obj.isMoving);
            initMoving = obj.isMoving; % store init
            for idx = movingIdx
                preWindow = initMoving(max(1, (idx - searchWindow)):idx-1);
                postWindow = initMoving((idx+1):min(length(initMoving), idx + searchWindow));
                obj.isMoving(idx) = ~((sum(preWindow) < peakWidth && sum(postWindow) < peakWidth));
            end
        end
            
        function out = bin_DFF(obj, data, heading)
            if nargin < 3
                data = obj.initializeData(obj.data, 'spikes');
                heading = obj.initializeData(obj.floating, 'heading');
            end
            binEdges = -180:obj.binWidth:180; % Because alpha is [-180,180];
            groups = discretize(heading, binEdges);
            u_groups = unique(groups);
            out = zeros(size(data,1), length(u_groups));
            for bin = 1:length(u_groups)
                for c = 1:size(data, 1)
                    out(c, bin) = mean(data(c, groups == u_groups(bin)));
                end
            end
        end
        
        function h = calculate_heading(obj,r,phi,alpha)
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

