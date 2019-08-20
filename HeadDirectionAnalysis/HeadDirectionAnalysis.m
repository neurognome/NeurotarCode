classdef HeadDirectionAnalysis < handle
    
    % Make sure you have the circular statistics toolbox
    properties (Constant = true)
        radius    double = 120; % Radius of the arena
    end
    
    properties
        data struct
        floating struct
        meanQuadrantCorrelation double
        directionInfo double
        
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
            
            obj.getHeading();
        end
        
        function out = getPlottingData(obj)
            out = obj.binData();  % Call w/o input arguments so it correctly gets the thing
        end
        
        function getPreferredDirection(obj, method)
            if nargin < 2
                method = 'vectorsum';
                fprintf('No method provided, defaulting to vectorsum\n')
            end
            
            binnedData = obj.binData();
            obj.directionInfo = zeros(size(binnedData, 1), 2);
            switch method
                case 'vectorsum'
                    % Quick anonymous functions
                    getHorz = @(v, theta) v .* cos(theta);
                    getVert = @(v, theta) v .* sin(theta);
                    getAng = @(vert, horz) atan(vert ./ horz);
                    getMag = @(vert, horz) sqrt(horz ^ 2 + vert ^ 2);
                    
                    for c = 1:size(binnedData, 1)
                        dat = binnedData(c, :);
                        dat(dat < 0) = 0; % rectify
                        theta = linspace(0, 2*pi, length(dat));
                        
                        % testing data : dat = circshift([9 8 7 6 5 6 7 8],ii);
                        h = getHorz(dat, theta);
                        v = getVert(dat, theta);
                        
                        r_h = sum(h);
                        r_v = sum(v);
                        
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
                        obj.directionInfo(c, :) = [ang, m];
                    end
                case 'max'
                    for c = 1:size(binnedData, 1)
                        [obj.directionInfo(c, 2) ,obj.directionInfo(c, 1)] = max(binnedData(c, :));
                    end
                otherwise
                    fprintf('Invalid method provided (methods: vectorsum or max)\n')
                    return
            end
        end    
        
        function calculateHeadDirectionIdx(obj, numSections)
            % the following analyses are based on Giacomo et al 2017, in
            % current biology...
            if nargin < 2
                numSections = 4;
            end
            
            neuralData = obj.initializeData(obj.data,'spikes');
            heading    = obj.initializeData(obj.floating,'heading');
            
            totalDur = length(heading);
            sectionLength = round(totalDur./numSections);
            
            q_data = cell(numSections, 1);
            q_heading = cell(numSections, 1);
            for q = 1:numSections % Separate timeseries into sections
                q_data{q}    = neuralData(:, ...
                    (q - 1) * sectionLength + 1 : min(q*sectionLength, length(neuralData)));
                q_heading{q} = heading(...
                    (q - 1) * sectionLength + 1 : min(q*sectionLength, length(heading)));
            end
            q_binnedData = zeros(size(neuralData, 1), length(-180:obj.binWidth:180)-1, numSections);
            % Calculate the binned data
            for q = 1:numSections
                q_binnedData(:, :, q) = obj.binData(q_data{q}, q_heading{q});
            end
            % compare quadrants in pairwise (1-1,1-2,1-3,1-4,2-3,2-4,3-4)
            possibleCombinations = nchoosek(1:4, 2);
            quadrantCorrelations = zeros(size(neuralData, 1), size(possibleCombinations, 1));
            for ii = 1:size(possibleCombinations, 1)
                for c = 1:size(neuralData, 1)
                    quadrantCorrelations(c, ii) = corr(...
                        q_binnedData(c, :, possibleCombinations(ii, 1))',...
                        q_binnedData(c, :, possibleCombinations(ii, 2))');
                end
            end
            obj.meanQuadrantCorrelation = mean(quadrantCorrelations,2);
        end
        
        function removeMovingSamples(obj, speedThreshold, peakWidth, searchWindow)
            if nargin < 2
                speedThreshold = 10; % Default
            end
            if nargin < 3
                peakWidth = 5;
            end
            if nargin < 4
                searchWindow = 10;
            end
            
            obj.findMovingTimes(speedThreshold, peakWidth, searchWindow);
            obj.removeSlowFlag = true;
        end
        
        %% Setters and Getters
        function setHeadingFlag(obj, val)
            obj.fixHeadingFlag = val;
            obj.getHeading();
        end
    end
    
    methods (Access = protected)
        function getHeading(obj)
            if obj.fixHeadingFlag
                obj.floating.heading = obj.calculateHeading(...
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
        
        function findMovingTimes(obj, speed_threshold, peakWidth, searchWindow)
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
        
        function out = binData(obj, data, heading)
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
        
        function h = calculateHeading(obj, r, phi,alpha)
            h = 2 .* asind((r .* sind(alpha - phi)) ./ (2 .* obj.radius)) + alpha;
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

