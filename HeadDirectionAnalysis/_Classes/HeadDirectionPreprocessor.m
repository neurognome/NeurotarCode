classdef HeadDirectionPreprocessor < NeurotarPreProcessor
    properties (Access = protected)
        flicker_flag = false;
    end
    
    methods
        function obj = HeadDirectionPreprocessor(data, floating, flicker_flag, varargin)
            obj@NeurotarPreProcessor(data, floating, varargin{:});
            obj.flicker_flag = flicker_flag;
            obj.checkForSpikes();
        end
        
        function checkForSpikes(obj)
            if ~isfield(obj.data, 'spikes')
                disp('Getting spikes')
                obj.data = obj.getSpikeData(obj.data);
            end
            obj.saveData();
        end
        
        function separateLightDarkData(obj, first_segment)
            if nargin < 2 || isempty(first_segment)
                first_segment = questdlg('Which segment came first?', 'First segment', 'light', 'dark', 'light');
            end
           
            segment_cell = inputdlg('How long was each segment (s):');
            
            
            segment_length = str2double(segment_cell{1}) * 10; % in frames
            recording_length = size(obj.data.DFF, 2);
            
            switch first_segment
                case 'light'
                    [light_data, dark_data] =  obj.separateRecording(obj.data, segment_length, recording_length);
                    [light_floating, dark_floating] =  obj.separateRecording(obj.floating, segment_length, recording_length);
                case 'dark'
                    [dark_data, light_data] = obj.separateRecording(obj.data, segment_length, recording_length);
                    [dark_floating, light_floating] =  obj.separateRecording(obj.floating, segment_length, recording_length);
            end
            
            mkdir('Separated Data')
            cd('Separated Data');
            
            save light_data.mat light_data
            save dark_data.mat dark_data
            save light_floating.mat light_floating
            save dark_floating.mat dark_floating
        end
        
        function [data1, data2] = separateRecording(obj, data, segment_length, recording_length)
            %% Get all the fields in the data structure
            
            data1 = struct();
            data2 = struct();
            struct_fields = fields(data);
            for i_field = 1:length(struct_fields)
                if ismember(recording_length, size(data.(struct_fields{i_field})))
                    [data1.(struct_fields{i_field}), data2.(struct_fields{i_field})] = ...
                        obj.separateField(data.(struct_fields{i_field}), segment_length, recording_length);
                else
                    data1.(struct_fields{i_field}) = data.(struct_fields{i_field});
                    data2.(struct_fields{i_field}) = data.(struct_fields{i_field});
                end
            end
        end
        
        function [data1, data2] = separateField(obj, data, segment_length, recording_length)
            % For separating out a single data
            
            n_segments = recording_length / segment_length;
            
            if mod(n_segments, 1) ~= 0 % Not an integer
                error('Something''s wrong, noninteger # of segments')
            end
            
            flip_flop = 1;
            data1 = [];
            data2 = [];
            for i_segment = 1:n_segments
                curr = (i_segment - 1) * segment_length;
                % Prepare data, this is finding the correct index to segment
                inds = repmat({1}, 1, ndims(data));
                for i_dim = 1:ndims(data)
                    inds{i_dim} = 1:size(data, i_dim);
                end
                
                working_dim = find(ismember(size(data), recording_length));
                inds{working_dim} = curr + 1: curr + segment_length;
                curr_data = data(inds{:});
                switch flip_flop
                    case 1
                        data1 = cat(working_dim, data1, curr_data);
                        flip_flop = 2;
                    case 2
                        data2 = cat(working_dim, data2, curr_data);
                        flip_flop = 1;
                end
            end
        end
        
        function processData(obj, averageWindow)
            processData@NeurotarPreProcessor(obj, averageWindow) % Run the superclass method
            if obj.flicker_flag
                obj.separateLightDarkData();
            end
        end
    end
    
    methods (Static)
        function data = getSpikeData(data)
            dffDeconv = zeros(size(data.DFF));
            for n = 1:size(data.DFF, 1)
                subroutine_progressbar(n/size(data.DFF,1));
                % get trace and run deconvolution
                trace = data.DFF(n,:);
                [~, spikes, ~] = deconvolveCa(trace, 'ar1' ,'foopsi', 'optimize_pars');
                dffDeconv(n, :) = spikes;
                
            end
            subroutine_progressbar(1);close all;
            
            data.spikes = dffDeconv;
        end
    end
    
end

