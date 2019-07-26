classdef NeurotarPreProcessor < handle
    
    % To preprocess neurotar data I guess...
    
    % looks like alpha is the angle the mouse is viewing, ie from the mouse head
    % omega is the posture angle, ie how much is the mouse tilted
    
    
    properties (Constant = true)
        variance_threshold(1,1) double = 0.001; % What threshold of variance within the neurotar sampling rate which is considered to be "good enough"
    end
    
    properties (Access = protected)
        imaging_data struct
        neurotar_data struct
    end
    
    properties (SetAccess = protected) % I have this as SetAccess protected because I want to see if it happened
        ForceTimeLock logical = 0
    end
    
    methods
        
        function obj = NeurotarPreProcessor(data,floating) % contstructor
            
            obj.imaging_data  = data;  % this lets our constructor assign the input arguments to the properties, for easier passing
            obj.neurotar_data = floating;
            
        end
        
        function processData(obj)

            time = obj.extractTime;
            
            if obj.ForceTimeLock
                isUniformSampling = false;
            else
                isUniformSampling = obj.checkSamplingRate(time);
            end
            
            obj.samplingFrequencyEqualizer(time,isUniformSampling);
            
            obj.preprocessChecker;
            
        end
        
        function obj = setForceTimeLock(obj, val)
            obj.ForceTimeLock = val;
        end
        
        function out = getImagingData(obj)
            out = obj.imaging_data;
        end
        
        function out = getNeurotarData(obj)
            out = obj.neurotar_data;
        end
        
    end
    
    
    methods (Access = protected)
        
        function time = extractTime(obj)
            time = obj.neurotar_data.time - '00:00:00.000'; % subtracting is necessary to turn the time into a matrix
            time = time(:, [4,5,7,8,10:12]); % this is assuming none of recording last more than an hour
            time = time .* [6 * 10^5, 6 * 10^4, 10^4, 10^3, 10^2, 10^1, 1];
            time = sum(time, 2) / 1000;
        end
        
        function isUniform = checkSamplingRate(obj,time)
            offset_times = cat(1,zeros(1,size(time,2)),time);
            deltaTime = time-offset_times(1:end-1,:);
            isUniform = var(deltaTime) < obj.variance_threshold;
        end
        
        function samplingFrequencyEqualizer(obj,time,isUniform)
            if isUniform
                fprintf('Your sampling rate is even, upsampling DFF to the neurotar frequency...\n')
                for ii =1:size(obj.imaging_data.DFF,1)
                    temp(ii,:) = resample(obj.imaging_data.DFF(ii,:),size(obj.neurotar_data.time,1)+1,size(obj.imaging_data.DFF,2));
                end
                obj.imaging_data.DFF = temp;
            else
                fprintf('Your sampling rate is too uneven, subsampling neurotar data to match DFF... \n')
                twoP_sampled_times = 1/obj.imaging_data.frameRate:1/obj.imaging_data.frameRate:(obj.imaging_data.numFrames / obj.imaging_data.frameRate);
                neurotar_matched_indices = zeros(1, length(twoP_sampled_times));
                for tt = 1:length(twoP_sampled_times)
                    [ ~, neurotar_matched_indices(tt) ] = min(abs(time - twoP_sampled_times(tt)));
                end
                obj.neurotar_data.X     = obj.neurotar_data.X(neurotar_matched_indices);
                obj.neurotar_data.Y     = obj.neurotar_data.Y(neurotar_matched_indices);
                obj.neurotar_data.phi   = obj.neurotar_data.phi(neurotar_matched_indices);
                obj.neurotar_data.alpha = obj.neurotar_data.alpha(neurotar_matched_indices);
                obj.neurotar_data.omega = obj.neurotar_data.omega(neurotar_matched_indices);
                obj.neurotar_data.speed = obj.neurotar_data.speed(neurotar_matched_indices);
            end
            
        end
        
        function preprocessChecker(obj)
            if length(obj.imaging_data.DFF) == length(obj.neurotar_data.X)
                fprintf('Your data are preprocessed, ready to go...\n')
            else
                fprintf('Something bad happened and your data are not yet compatible.\n')
                return
            end
        end
        
        
    end
    
end

