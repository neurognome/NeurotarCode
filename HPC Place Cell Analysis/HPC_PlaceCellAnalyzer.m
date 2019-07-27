classdef HPC_PlaceCellAnalyzer < handle
    
    
    properties (Access = public) %protected
        imaging_data struct
        neurotar_data struct
        template struct
        
        numBins double = 10
    end
    
    properties (Access = public ) %private
        workingData struct % this is an easy catchall struct that just puts all our data together for easier management
        DFF_binned double
        heatmaps double
        
        % none of these actually need to be properties... they're just used in interim for calculating things.
        
        spatial_info double
        isPlaceCell logical
    end
    
    methods
        function obj = HPC_PlaceCellAnalyzer(data,floating) % Contsructor funtion
            obj.imaging_data = data;
            obj.neurotar_data = floating;
            
            obj.importTemplate(); % Initial import of the template
            [obj.workingData] = binDFF(obj); % Initialize working data by passing in "working variables". These are variables which may be used or not...
        end
        
        
        function findPlaceCells(obj)
            obj.computeSpatialInformation()
            obj.screenPlaceCells()
            
            % plot
            figure
            N = length(obj.spatial_info);
            plot(1:N, obj.spatial_info, 'k-'); hold on
            plot(find(obj.isPlaceCell == 1), obj.spatial_info(find(obj.isPlaceCell == 1)), 'ro', 'MarkerFaceColor', 'r'); hold on
            plot([1, N], [0, 0], 'k--')
        end
        
        function makeHeatMaps(obj,cellNum) % i only need this to run on a single cell, which I've specfied
            obj.computeSpatialInformation()
            
            imagesc(obj.heatmaps(:,:,cellNum)); hold on
            title(strcat({'Heat map of neuron'}, {' '},  num2str(cellNum)));
            axis square
        end
        
        
        function computeSpatialInformation(obj)
            unpackData(obj.workingData,'counts')
            for ii = 1:size(obj.DFF_binned,3)
                curr_dat = obj.DFF_binned(:,:,ii);
                obj.heatmaps = counts / sum(counts(:)) .* curr_dat / nanmean(curr_dat(:)) .* log2(curr_dat / nanmean(curr_dat(:)));
                obj.heatmaps(isnan(obj.heatmaps)) = 0;
                obj.spatial_info(ii) = sum(obj.heatmaps(:));
            end
        end % lifted from Will's code, Step1/2
        
        function screenPlaceCells(obj) % lifted from Will's code Step2/2
            unpackData(obj.workingData,'counts','bin_id_X','bin_id_Y','cell_responses')
            
            for ii = 1:length(obj.spatial_info)
                if obj.spatial_info(ii) > 0
                    num_shuffles = 500;
                    shuff_spatial_info = zeros(1, num_shuffles);
                    for jj = 1:num_shuffles
                        shuff_DFF_binned = zeros(size(obj.DFF_binned,1),size(obj.DFF_binned,2));
                        shuff_start_point = randi(length(cell_responses(ii,:)));
                        shuff_cell_responses = cell_responses(ii,[shuff_start_point:length(cell_responses(ii,:)),...
                            1:(shuff_start_point - 1)]);
                        for kk = 1:size(bin_id_X,1)
                            shuff_DFF_binned(bin_id_X(kk,ii), bin_id_Y(kk,ii)) = shuff_DFF_binned(bin_id_X(kk,ii), bin_id_Y(kk,ii)) + ...
                                shuff_cell_responses(kk);
                        end
                        shuff_DFF_binned = shuff_DFF_binned ./ counts;
                        
                        I_shuff = counts / sum(counts(:)) .* shuff_DFF_binned / nanmean(shuff_DFF_binned(:)) .* ...
                            log2(shuff_DFF_binned / nanmean(shuff_DFF_binned(:)));
                        I_shuff(isnan(I_shuff)) = 0;
                        shuff_spatial_info(jj) = sum(I_shuff(:));
                    end
                    
                    percentile = length(find(obj.spatial_info(ii) > shuff_spatial_info)) / num_shuffles;
                    if percentile > 0.8
                        obj.isPlaceCell(ii) = 1;
                    else
                        obj.isPlaceCell(ii) = 0;
                    end
                else
                    obj.isPlaceCell(ii) = 0;
                end
            end
        end
        
        function out = binDFF(obj)
            bin_X = obj.template.X(1):(obj.template.X(end) - obj.template.X(1))/(obj.numBins - 1):obj.template.X(end);
            bin_Y = obj.template.Y(1):(obj.template.Y(end) - obj.template.Y(1))/(obj.numBins - 1):obj.template.Y(end);
            
            cell_responses = (obj.imaging_data.DFF - min(obj.imaging_data.DFF , [], 2)) ./ ...
                (max(obj.imaging_data.DFF , [], 2) - min(obj.imaging_data.DFF , [], 2)); % normalizing the cell responses to make the heat map more interpretible
            
            N = size(cell_responses,1);
            bin_id_X = zeros(length(obj.neurotar_data.X),N);
            bin_id_Y = zeros(length(obj.neurotar_data.Y),N);
            
            for ii = 1:N
                disp(num2str(ii))
                dff_temp = zeros(obj.numBins, obj.numBins);
                ct       = zeros(obj.numBins,obj.numBins);
                
                for jj = 1:length(obj.neurotar_data.X)
                    x_distance = abs(obj.neurotar_data.X(jj) - bin_X);
                    y_distance = abs(obj.neurotar_data.Y(jj) - bin_Y);
                    [~, bin_id_X(jj,ii)] = min(x_distance);
                    [~, bin_id_Y(jj,ii)] = min(y_distance);
                    
                    ct(bin_id_X(jj,ii), bin_id_Y(jj,ii)) = ct(bin_id_X(jj,ii), bin_id_Y(jj,ii)) + 1;
                    
                    dff_temp(bin_id_X(jj,ii), bin_id_Y(jj,ii)) = dff_temp(bin_id_X(jj,ii), bin_id_Y(jj,ii)) + cell_responses(ii, jj);
                end
                obj.DFF_binned(:,:,ii) = dff_temp;
            end
            
            counts = ct; %doesn't change
            obj.DFF_binned = obj.DFF_binned ./ counts;
            
            out = packData(counts,bin_id_X,bin_id_Y,cell_responses);
        end % this code is lifted from Will's code, just some minor changes to make it compatible
        
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
end

