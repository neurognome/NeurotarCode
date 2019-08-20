classdef NeurotarPreProcessor < handle
    
    % To preprocess neurotar data I guess...
    
    % looks like alpha is the angle the mouse is viewing, ie from the mouse head
    % omega is the posture angle, ie how much is the mouse tilted
    
    
    properties (Constant = true)
        variance_threshold double = 0.001; % What threshold of variance within the neurotar sampling rate which is considered to be "good enough"
    end
    
    properties (Access = protected)
        data 
        floating
        
        initData
        initFloating
    end
    
    properties (SetAccess = protected) % I have this as SetAccess protected because I want to see if it happened
        forceTimeLock logical = 0
        averageWindow double = 0; % Window for averaging for time lock stimulus
   end
    
    methods
        
        function obj = NeurotarPreProcessor(data,floating) % contstructor
            
            if nargin == 0
                fprintf('Choose your data file: \n')
                data = importdata(uigetfile('.mat'));
                fprintf('Choose your stimulus file: \n')
                floating = importdata(uigetfile('.mat'));
            end
            obj.data = data;
            obj.floating = floating;
            
            obj.initData  = data;  % This stores the initial data away in case we need to revert to it
            obj.initFloating = floating;
        end
        
        function [data, floating] = processData(obj, averageWindow)    
            if exist('averageWindow', 'var')  % If supplied average window, else default
               obj.setAverageWindow(averageWindow);
            end 
            
            time = obj.extractTime(obj.floating);  % Extract times from neurotar
            
            if obj.forceTimeLock  % Upsample DFF or downsample neurotar? This forces downsampling
                isUniformSampling = false;
            else
                isUniformSampling = obj.checkSamplingRate(time);
            end
            
            obj.equalizeFs(time, isUniformSampling);
            
            obj.preprocessChecker();  % Confirming everything worked by checking lengths against each other
            
            data = obj.getData();
            floating = obj.getFloating();
        end
        
        function reset(obj)
            obj.data = obj.initData;
            obj.floating = obj.initFloating;
            fprintf('Data reset to initial values/n')
        end
        
        %% Getters n Setters
        function setAverageWindow(obj, win)
            if mod(win, 2) == 1
                fprintf('Window must be even, adjusted: %d --> %d\n', win, win - 1);
                win = win - 1;
            end
            obj.averageWindow = win;
        end
        
        function obj = setForceTimeLock(obj, val)
            obj.forceTimeLock = val;
        end
        
        function out = getData(obj)
            out = obj.data;
        end
        
        function out = getFloating(obj)
            out = obj.floating;
        end
    end
    
    
    methods (Access = protected)
        function time = extractTime(obj, floating)
            time = floating.time - '00:00:00.000'; % subtracting is necessary to turn the time into a matrix
            time = time(:, [4,5,7,8,10:12]); % this is assuming none of recording last more than an hour
            time = time .* [6 * 10^5, 6 * 10^4, 10^4, 10^3, 10^2, 10^1, 1];
            time = sum(time, 2) / 1000;
        end
        
        function isUniform = checkSamplingRate(obj,time)
            offset_times = cat(1,zeros(1,size(time,2)),time);
            deltaTime = time-offset_times(1:end-1,:);
            isUniform = var(deltaTime) < obj.variance_threshold;
        end
        
        function equalizeFs(obj, time, isUniform)
            if isUniform
                fprintf('Your sampling rate is even, upsampling DFF to the neurotar frequency...\n')
                
                for ii =1:size(obj.data.DFF,1)
                    obj.data.DFF(ii,:) = resample(obj.data.DFF(ii,:),size(obj.floating.time,1)+1,size(obj.data.DFF,2));
                end
                                
                if isfield(obj.data,'spikes')
                    for ii =1:size(obj.data.spikes,1)
                        obj.data.spikes(ii,:) = resample(obj.data.spikes(ii,:),size(obj.floating.time,1)+1,size(obj.data.spikes,2));
                    end
                end
                
            else
                if obj.forceTimeLock
                    fprintf('Forced to lock to stimulus time, not recommended... \n')
                else
                    fprintf('Your sampling rate is too uneven, subsampling neurotar data to match DFF... \n')
                end
                twoP_sampled_times = 1/obj.data.frameRate:1/obj.data.frameRate:(obj.data.numFrames / obj.data.frameRate);
                neurotar_matched_indices = zeros(1, length(twoP_sampled_times));
                for tt = 1:length(twoP_sampled_times)
                    [ ~, neurotar_matched_indices(tt) ] = min(abs(time - twoP_sampled_times(tt)));
                    
                   %timeDiff = round(abs(time - twoP_sampled_times(tt)), 5);
                  % neurotar_matched_indices(tt) = min(find(timeDiff == min(timeDiff))); % These lines of code are just 
                   % so that we can match python. See the python version for details
                end
                obj.extractFloatingFields(neurotar_matched_indices);
            end
 
        end
        
        function extractFloatingFields(obj, idx)
            floatingFields = string(fields(obj.floating))';  % Convert into a "string array"
            for f = floatingFields
                if strcmp(f, 'time')
                    continue
                else
                    initField = obj.floating.(f);  % This is because if we change floating iteratively, it'll affect the 
                    % next window as well
                    newField = zeros(length(idx), 1);
                    for ii = 1:length(idx)
                       newField(ii) = ...
                            mean(initField(...
                            max([1, idx(ii) - floor(obj.averageWindow / 2)]) : ...
                            min([idx(ii) + floor(obj.averageWindow / 2), length(obj.floating.(f))])));
                    end
                    obj.floating.(f) = newField;
                end
            end
        end
        
        function preprocessChecker(obj)
            if length(obj.data.DFF) == length(obj.floating.X)
                fprintf('Your data are preprocessed, ready to go...\n')
            else
                fprintf('Something bad happened and your data are not yet compatible.\n')
                return
            end
        end
        
        
    end
    
end

