classdef NeurotarDataExtractor < handle
    %-------------------------------------------------------------------------%
    %   This script reads the excel file generated by the Neurotar (make sure
    %   you've downloaded the excel converter), saves the relevant spatial
    %   parameters, and then computes a occupancy map and percentage of time
    %   the mouse spent moving.
    %
    %   Written by WTR 05/10/2019 // Last updated by WTR 07/18/2019
    %
    %   Converted into class: KS 31Jul2019
    %-------------------------------------------------------------------------%
    
    properties
        mouse string
        date  string
        image_flag logical
        session_num double
        extractedData
    end
    
    methods
        function obj = NeurotarDataExtractor(mouse,date,session_num,image_flag)
            if nargin < 4
                args = obj.getUserInput();
                obj.mouse = args{1};
                obj.date  = args{2};
                obj.session_num = str2double(args{3});
                obj.image_flag   = logical((args{4})); % because it doesn't come in as a bool
            else
                obj.mouse = mouse;
                obj.date = date;
                obj.session_num = session_num;
                obj.image_flag = image_flag;
            end
            
            if obj.image_flag
                filename = strcat(obj.mouse, '_', obj.date, '_IMAGING_', num2str(obj.session_num), '.xlsx');
            else
                filename = strcat(obj.mouse, '_', obj.date, '.xlsx');
            end
            [~, ~, behaviour_data] = xlsread(filename, 3);
            
            % Initial methods
            obj.extractData(behaviour_data);
        end
        
        function saveData(obj)
            floating = obj.extractedData.export();
            save(strcat('floating_data_', obj.mouse, '_', obj.date, '_session', num2str(obj.session_num), '.mat'), 'floating');
        end
        
        
        function visualize(obj)
            obj.extractedData.exportVar(); % Dirty way of getting compatibility
            
            %% Occupancy distribution
            % Pulling the relevant data structures from the Neurotar excel sheet.
            
            max_X = max(X); min_X = min(X);
            max_Y = max(Y); min_Y = min(Y);
            
            diff_time = time(2:end, :) - time(1:(end - 1), :);
            diff_time = diff_time(:, [4,5,7,8,10:12]); % this is assuming none of recording last more than an hour
            diff_time_in_millisecs = diff_time .* [60 * 10^5, 60 * 10^4, 10^4, 10^3, 10^2, 10^1, 1];
            delta_T = sum(diff_time_in_millisecs, 2) / 1000; % delta_T is in seconds
            
            % Removing artifacts from when the magnets are off the table and the system has no idea where the mouse is
            artifacts = find(speed > 300);
            speed(artifacts) = [];
            
            % Computing the occupancy map and correlation between Neurotar's outputted
            % speed and a - non-filtered - speed from the position alone.
            n_bins = 100;
            bin_X = linspace(min_X, max_X, n_bins);
            bin_Y = linspace(min_Y, max_Y, n_bins);
            counts = zeros(n_bins, n_bins);
            zoned_bins = zeros(n_bins, n_bins);
            inst_speed = zeros(1, length(X) - 1);
            
            for ii = 2:(length(X) - 2)
                x_distance = abs(X(ii) - bin_X);
                y_distance = abs(Y(ii) - bin_Y);
                [~, bin_id_X] = min(x_distance);
                [~, bin_id_Y] = min(y_distance);
                counts(bin_id_X, bin_id_Y) = counts(bin_id_X, bin_id_Y) + 1;
                zoned_bins(bin_id_X, bin_id_Y) = zoned_bins(bin_id_X, bin_id_Y) + zones(ii);
                
                dX = X(ii + 1) - X(ii);
                dY = Y(ii + 1) - Y(ii);
                inst_speed(ii - 1) = sqrt(dX^2 + dY^2) / delta_T(ii);
                
            end
            
            inst_speed(artifacts) = [];
            
            % Plotting
            figure
            image(counts);
            axis square
            title('Distribution of visits');
            
            zoned_bins = zoned_bins ./ counts;
            zoned_bins(isnan(zoned_bins)) = 0;
            
            figure
            imagesc(zoned_bins);
            axis square
            title('Zone map');
            
            %% Speed distribution
            figure
            histogram(speed);
            title('Speed distribution');
            xlabel('Speed (mm^2/sec)');
            
            figure
            plot(1:length(speed), speed, 'k-'); hold on
            title('Speed trace');
            xlabel('Sampling times');
            ylabel('Speed (mm^2 / sec');
            
            % speed_corr = corr(speed(2:end), inst_speed')
            norm_speed = (speed - min(speed)) / (max(speed) - min(speed)); % @ Joe's suggestion
            time_moving = length(norm_speed(norm_speed > 0.10)) / length(norm_speed)
        end
    end
    
    methods (Access = protected)
        function out = getUserInput(obj)
            
            prompt = {'Enter mouse name: ',...
                'Enter recording date: ',...
                'Enter session number: ',...
                'Imaged? (yes = 1)'};
            
            dlgtitle = 'User input';
            
            definput = obj.getDefaultInput;
            
            out = inputdlg(prompt,dlgtitle,[1 35],definput);
            
        end
        
        function out = getDefaultInput(obj)
            file = dir('*.xlsx');
            if length(file) > 1
                disp('Choose your excel file...')
                file = dir(uigetfile('*.xlsx')); % Get the right one
            end
            
            mouse       = file.name(1:6);
            date        = file.name(8:12);
            session_num = file.name(22);
            if ~isempty(strfind(file.name,'IMAGING'))
                imaged = '1';
            else
                imaged = '0';
                
            end
            out = {mouse, date, session_num, imaged};
        end
        
        function extractData(obj,behaviour_data)
            obj.extractedData = DataObject();
            X     = cell2mat(behaviour_data(2:end, 7)); % their X and Y are flipped on their data sheet
            Y     = cell2mat(behaviour_data(2:end, 6));
            r     = cell2mat(behaviour_data(2:end,3));
            phi   = cell2mat(behaviour_data(2:end, 4));
            alpha = cell2mat(behaviour_data(2:end, 5));
            omega = cell2mat(behaviour_data(2:end, 8));
            speed = cell2mat(behaviour_data(2:end, 9));
            time  = cell2mat(behaviour_data(3:end, 2));
            zones = cell2mat(behaviour_data(2:end, 10));
            
            obj.extractedData.add('X','Y','r','phi','alpha','omega','speed','time','zones');
        end
        
    end
    
end



