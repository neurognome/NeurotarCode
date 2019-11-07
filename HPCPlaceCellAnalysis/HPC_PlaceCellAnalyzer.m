classdef HPC_PlaceCellAnalyzer < handle
    % ------------------------------------------------------------------------
    % For analyzing hippocampal place cell data from Neurotar experiments.
    % 90% of the code is taken from Will's code, and just reformatted to work
    % in this configuration.
    %
    % Dependencies:
    % dataObject <https://github.com/kevinksit/GeneralHelperCode>
    %
    % Written 27Jul2019 KS
    % Updated
    % ------------------------------------------------------------------------
    properties (Constant = true)
        numBins double = 10
    end
    
    properties (SetAccess = protected)
        imaging_data struct
        neurotar_data struct
        template struct
        
        % Hierarchy
        analysisData
        plottingData
        workingData
        
    end
    
    %% these are public methods, which allows us to call them from a script or the command window
    methods
        function obj = HPC_PlaceCellAnalyzer(data,floating) % Contsructor funtion
            obj.imaging_data = data;
            obj.neurotar_data = floating;
            
            % Initialize data objects
            obj.workingData = DataObject();
            obj.plottingData = DataObject();
            obj.analysisData = DataObject();
            
            obj.importTemplate(); % Initial import of the template
            obj.binDFF(); % Initialize working data by passing in "working variables". These are variables which may be used or not...
        end
        
        function findPlaceCells(obj)
            obj.computeSpatialInformation()
            obj.screenPlaceCells()
            
            % plot
            figure
            N = length(obj.analysisData.spatial_info);
            plot(1:N, obj.analysisData.spatial_info, 'k-'); hold on
            plot(find(obj.analysisData.isPlaceCell == 1), obj.analysisData.spatial_info(find(obj.analysisData.isPlaceCell == 1)), 'ro', 'MarkerFaceColor', 'r'); hold on
            plot([1, N], [0, 0], 'k--')
        end
        
        function makeHeatMap(obj,cellNum) % i only need this to run on a single cell, which I've specfied
            obj.computeSpatialInformation()
            
            figure
            imagesc(obj.plottingData.I(:,:,cellNum)); hold on
            title(strcat({'Heat map of neuron'}, {' '},  num2str(cellNum)));
            axis square
        end
        
        function changeTemplate(obj) % This lets you change your template, then re-bin according to the new templaet..
            fprintf('Choose your new template...\n')
            [fn,pn] = uigetfile('.mat');
            obj.importTemplate(pn,fn);
        end
    end
    
    
    %% The following methods are protected, ie they can't be accessed outside of this object. They're for internal use only
    
    methods (Access = protected)
        function binDFF(obj)
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
                DFF_binned(:,:,ii) = dff_temp;
            end
            
            counts = ct; %doesn't change
            DFF_binned = DFF_binned ./ counts;
            
            obj.analysisData.add('DFF_binned');
            obj.workingData.add('counts','bin_id_X','bin_id_Y','cell_responses');
        end % this code is lifted from Will's code, just some minor changes to make it compatible
        
        function computeSpatialInformation(obj)
            for ii = 1:size(obj.analysisData.DFF_binned,3)
                curr_dat = obj.analysisData.DFF_binned(:,:,ii);
                heatmaps = obj.workingData.counts / sum(obj.workingData.counts(:)) .* curr_dat / nanmean(curr_dat(:)) .* log2(curr_dat / nanmean(curr_dat(:)));
                heatmaps(isnan(heatmaps)) = 0;
                spatial_info(ii) = sum(heatmaps(:));
                I(:,:,ii) = heatmaps;
            end
            
            obj.analysisData.add('spatial_info');
            obj.plottingData.add('I');
        end % lifted from Will's code, Step1/2
        
        function screenPlaceCells(obj) % lifted from Will's code Step2/2
            % Although not totally recommended, I exported these two because they're accessed so much in the shuffling
            % that keeping them in the object was slowing down the processing a LOT
            obj.workingData.exportData('bin_id_X')
            obj.workingData.exportData('bin_id_Y')
            
            for ii = 1:length(obj.analysisData.spatial_info)
                if obj.analysisData.spatial_info(ii) > 0
                    num_shuffles = 500;
                    shuff_spatial_info = zeros(1, num_shuffles);
                    for jj = 1:num_shuffles
                        shuff_DFF_binned = zeros(size(obj.analysisData.DFF_binned,1),size(obj.analysisData.DFF_binned,2));
                        shuff_start_point = randi(length(obj.workingData.cell_responses(ii,:)));
                        shuff_cell_responses = obj.workingData.cell_responses(ii,[shuff_start_point:length(obj.workingData.cell_responses(ii,:)),...
                            1:(shuff_start_point - 1)]);
                        for kk = 1:size(bin_id_X,1)
                            shuff_DFF_binned(bin_id_X(kk,ii), bin_id_Y(kk,ii)) = shuff_DFF_binned(bin_id_X(kk,ii), bin_id_Y(kk,ii)) + ...
                                shuff_cell_responses(kk);
                        end
                        shuff_DFF_binned = shuff_DFF_binned ./ obj.workingData.counts;
                        
                        I_shuff = obj.workingData.counts / sum(obj.workingData.counts(:)) .* shuff_DFF_binned / nanmean(shuff_DFF_binned(:)) .* ...
                            log2(shuff_DFF_binned / nanmean(shuff_DFF_binned(:)));
                        I_shuff(isnan(I_shuff)) = 0;
                        shuff_spatial_info(jj) = sum(I_shuff(:));
                    end
                    
                    percentile = length(find(obj.analysisData.spatial_info(ii) > shuff_spatial_info)) / num_shuffles;
                    if percentile > 0.8
                        isPlaceCell(ii) = 1;
                    else
                        isPlaceCell(ii) = 0;
                    end
                else
                    isPlaceCell(ii) = 0;
                end
            end
            obj.analysisData.addData('isPlaceCell');
        end
                
        function importTemplate(obj,pn,fn)
            if nargin < 2
                fprintf('No template provided, loading default empty template... \n')
                obj.template = importdata('D:\Dropbox\Floating Cage\Arena templates\FloatCage_Empty_TEMPLATE.mat');
            else
                obj.template = load([pn fn]);
            end
        end
        
    end
end

