classdef HPC_PlaceCellAnalyzer < handle
    
    
    properties (Access = public) %protected
        imaging_data struct
        neurotar_data struct
        template struct
        
        numBins double = 10
    end
    
    properties (Access = public ) %private
        counts double
        DFF_binned double
        cell_responses double
        bin_id_X double
        bin_id_Y double
        
        spatial_info double

        isPlaceCell logical
    end
    
    methods
        function obj = HPC_PlaceCellAnalyzer(data,floating) % Contsructor funtion
            obj.imaging_data = data;
            obj.neurotar_data = floating;
            
            obj.importTemplate; % Initial import of the template
            
            obj.binDFF;
        end
        
        
        function out = findPlaceCells(obj);
            obj.computeSpatialInformation
            out = obj.screenPlaceCells
            
            %plot the stuff
        end
        
        function makeHeatMaps(obj,cellNum) % i only need this to run on a single cell, which I've specfied
            obj.makeBins;
            obj.binDFF;
            obj.computeSpatialInformation
            
            imagesc the hetamaps (I)
        end

        
        function computeSpatialInformation(obj)
            for ii = 1:size(obj.DFF_binned,3)
                curr_dat = obj.DFF_binned(:,:,ii);
                I = obj.counts / sum(obj.counts(:)) .* curr_dat / nanmean(curr_dat(:)) .* log2(curr_dat / nanmean(curr_dat(:)));
                I(isnan(I)) = 0;
                obj.spatial_info(ii) = sum(I(:));
            end
        end
        
        function [place_cell] = screenPlaceCells(obj)
            for ii = 1:length(obj.spatial_info)
                if obj.spatial_info(ii) > 0
                    num_shuffles = 500;
                    shuff_spatial_info = zeros(1, num_shuffles);
                    for jj = 1:num_shuffles
                        shuff_DFF_binned = zeros(size(obj.DFF_binned,1),size(obj.DFF_binned,2));
                        shuff_start_point = randi(length(obj.cell_responses(ii,:)));
                        shuff_cell_responses = obj.cell_responses(ii,[shuff_start_point:length(obj.cell_responses(ii,:)),...
                            1:(shuff_start_point - 1)]);
                        for kk = 1:size(obj.bin_id_X,1)
                            shuff_DFF_binned(obj.bin_id_X(kk,ii), obj.bin_id_Y(kk,ii)) = shuff_DFF_binned(obj.bin_id_X(kk,ii), obj.bin_id_Y(kk,ii)) + ...
                                shuff_cell_responses(kk);
                        end
                        shuff_DFF_binned = shuff_DFF_binned ./ obj.counts;
                        
                        I_shuff = obj.counts / sum(obj.counts(:)) .* shuff_DFF_binned / nanmean(shuff_DFF_binned(:)) .* ...
                            log2(shuff_DFF_binned / nanmean(shuff_DFF_binned(:)));
                        I_shuff(isnan(I_shuff)) = 0;
                        shuff_spatial_info(jj) = sum(I_shuff(:));
                    end
                    
                    percentile = length(find(obj.spatial_info(ii) > shuff_spatial_info)) / num_shuffles;
                    disp(num2str(percentile))
                    if percentile > 0.8
                        isPlaceCell(ii) = 1;
                    else
                        isPlaceCell(ii) = 0;
                    end
                else
                    isPlaceCell(ii) = 0;
                end
            end
        end
        
        function plotPlaceCells(obj)
        figure
        plot(1:N, obj.spatial_info, 'k-'); hold on
        plot(find(isPlaceCell == 1), obj.spatial_info(find(isPlaceCell == 1)), 'ro', 'MarkerFaceColor', 'r'); hold on
        plot([1, N], [0, 0], 'k--')
        end


        function binDFF(obj)
            bin_X = obj.template.X(1):(obj.template.X(end) - obj.template.X(1))/(obj.numBins - 1):obj.template.X(end);
                   bin_Y = obj.template.Y(1):(obj.template.Y(end) - obj.template.Y(1))/(obj.numBins - 1):obj.template.Y(end);
            
            obj.cell_responses = (obj.imaging_data.DFF - min(obj.imaging_data.DFF , [], 2)) ./ ...
                (max(obj.imaging_data.DFF , [], 2) - min(obj.imaging_data.DFF , [], 2)); % normalizing the cell responses to make the heat map more interpretible
            
            N = size(obj.cell_responses,1);
            obj.bin_id_X = zeros(length(obj.neurotar_data.X),N);
            obj.bin_id_Y = zeros(length(obj.neurotar_data.Y),N);
            
            for ii = 1:N
                disp(num2str(ii))
                dff_temp = zeros(obj.numBins, obj.numBins);
                bin_id_X = zeros(1, length(obj.neurotar_data.X));
                bin_id_Y = zeros(1, length(obj.neurotar_data.Y));
                ct       = zeros(obj.numBins,obj.numBins);
                
                for jj = 1:length(obj.neurotar_data.X)
                    x_distance = abs(obj.neurotar_data.X(jj) - bin_X);
                    y_distance = abs(obj.neurotar_data.Y(jj) - bin_Y);
                    [~, bin_id_X(jj)] = min(x_distance);
                    [~, bin_id_Y(jj)] = min(y_distance);
                    
                    ct(bin_id_X(jj), bin_id_Y(jj)) = ct(bin_id_X(jj), bin_id_Y(jj)) + 1;
                    
                    dff_temp(bin_id_X(jj), bin_id_Y(jj)) = dff_temp(bin_id_X(jj), bin_id_Y(jj)) + obj.cell_responses(ii, jj);
                end
                obj.DFF_binned(:,:,ii) = dff_temp;
                obj.bin_id_X(:,ii) = squeeze(bin_id_X);
                obj.bin_id_Y(:,ii) = squeeze(bin_id_Y);
            end
            
            obj.counts = ct; %doesn't change
            obj.DFF_binned = obj.DFF_binned ./ obj.counts;
        end
        
        
    end

     function importTemplate(obj,pn,fn)
            if nargin < 2
                fprintf('No template provided, loading default empty template... \n')
                obj.template = importdata('D:\Dropbox\Floating Cage\Arena templates\FloatCage_Empty_TEMPLATE.mat');
            else
                obj.template = load([pn fn]);
            end
        end
        
         function changeTemplate(obj) % This lets you change your template, then re-bin according to the new templaet..
            fprintf('Choose your new template...\n')
            [fn,pn] = uigetfile('.mat');
            obj.importTemplate(pn,fn);
        end
    
end


%{


            %% Step 1: Computing the spatial information
            I = counts / sum(counts(:)) .* DFF_binned / nanmean(DFF_binned(:)) .* log2(DFF_binned / nanmean(DFF_binned(:)));
            I(isnan(I)) = 0;
            spatial_info = sum(I(:));
            
            %% Step 2: Shuffling
            if spatial_info > 0
                num_shuffles = 500;
                shuff_spatial_info = zeros(1, num_shuffles);
                for jj = 1:num_shuffles
                    shuff_DFF_binned = zeros(size(DFF_binned));
                    shuff_start_point = randi(length(cell_responses));
                    shuff_cell_responses = cell_responses([shuff_start_point:length(cell_responses),...
                        1:(shuff_start_point - 1)]);
                    for kk = 1:length(bin_X)
                        shuff_DFF_binned(bin_X(kk), obj.bin_Y(kk)) = shuff_DFF_binned(bin_X(kk), bin_Y(kk)) + ...
                            shuff_cell_responses(kk);
                    end
                    shuff_DFF_binned = shuff_DFF_binned ./ counts;
                    
                    I_shuff = counts / sum(counts(:)) .* shuff_DFF_binned / nanmean(shuff_DFF_binned(:)) .* ...
                        log2(shuff_DFF_binned / nanmean(shuff_DFF_binned(:)));
                    I_shuff(isnan(I_shuff)) = 0;
                    shuff_spatial_info(jj) = sum(I_shuff(:));
                end
                
                percentile = length(find(spatial_info > shuff_spatial_info)) / num_shuffles;
                
                if percentile > 0.8
                    place_cell = 1;
                else
                    place_cell = 0;
                end
            else
                place_cell = 0;
            end
        end
        
    end
    end
%}