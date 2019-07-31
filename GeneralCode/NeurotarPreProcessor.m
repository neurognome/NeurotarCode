classdef NeurotarPreProcessor < handle
    
    % To preprocess neurotar data I guess...
    
    % looks like alpha is the angle the mouse is viewing, ie from the mouse head
    % omega is the posture angle, ie how much is the mouse tilted
    
    
    properties (Constant = true)
        variance_threshold(1,1) double = 0.001; % What threshold of variance within the neurotar sampling rate which is considered to be "good enough"
    end
    
    properties (Access = protected)
        initData
    end
    
    properties (SetAccess = protected) % I have this as SetAccess protected because I want to see if it happened
        ForceTimeLock logical = 0
        workingData
        
    end
    
    methods
        
        function obj = NeurotarPreProcessor(data,floating) % contstructor
            obj.workingData = dataObject();
            obj.initData  = dataObject(data,floating);  % This stores the initial data away in case we need to revert to it
            
            obj.processData(data,floating);
        end
        
        
        function processData(obj,data,floating)
            if nargin < 2
                data = obj.initData.data;
                floating = obj.initData.floating;
                
                obj.workingData.clearAllData();
            end
            
            time = obj.extractTime(floating);
            
            if obj.ForceTimeLock
                isUniformSampling = false;
            else
                isUniformSampling = obj.checkSamplingRate(time);
            end
            
            [data,floating] = obj.samplingFrequencyEqualizer(data,floating,time,isUniformSampling);
            
            obj.preprocessChecker(data,floating);
            
            obj.workingData = dataObject(data,floating);             % after all the processing steps, assign to working data. This can be called or further passed for more stuff if necessary..
        end
        
        function obj = setForceTimeLock(obj, val)
            obj.ForceTimeLock = val;
        end
    end
    
    
    methods (Access = protected)
        
        
        function time = extractTime(obj,floating)
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
        
        function [data,floating] = samplingFrequencyEqualizer(obj,data,floating,time,isUniform)
            if isUniform
                fprintf('Your sampling rate is even, upsampling DFF to the neurotar frequency...\n')
                for ii =1:size(data.DFF,1)
                    temp(ii,:) = resample(data.DFF(ii,:),size(floating.time,1)+1,size(data.DFF,2));
                end
                data.DFF = temp;
            else
                if obj.ForceTimeLock
                    fprintf('Forced to lock to stimulus time, not recommended... \n')
                else
                    fprintf('Your sampling rate is too uneven, subsampling neurotar data to match DFF... \n')
                end
                twoP_sampled_times = 1/data.frameRate:1/data.frameRate:(data.numFrames / data.frameRate);
                neurotar_matched_indices = zeros(1, length(twoP_sampled_times));
                for tt = 1:length(twoP_sampled_times)
                    [ ~, neurotar_matched_indices(tt) ] = min(abs(time - twoP_sampled_times(tt)));
                end
                
                floating.X     = floating.X(neurotar_matched_indices);
                floating.Y     = floating.Y(neurotar_matched_indices);
                floating.phi   = floating.phi(neurotar_matched_indices);
                floating.alpha = floating.alpha(neurotar_matched_indices);
                floating.omega = floating.omega(neurotar_matched_indices);
                floating.speed = floating.speed(neurotar_matched_indices);
            end
            
        end
        
        function preprocessChecker(obj,data,floating)
            if length(data.DFF) == length(floating.X)
                fprintf('Your data are preprocessed, ready to go...\n')
            else
                fprintf('Something bad happened and your data are not yet compatible.\n')
                return
            end
        end
        
        
    end
    
end

