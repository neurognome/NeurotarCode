classdef HPC_PlaceCellAnalyzer < handle
    

properties (Access = protected)
    imaging_data struct
    neurotar_data struct
    template struct

    numBins double = 10
end

properties (Access = private)
    counts double
    binned_DFF double
    bin_X double
    bin_Y double
    cell_responses double
    spatial_info double
end

methods
    function obj = HPC_PlaceCellAnalyzer(data,floating) % Contsructor funtion
        obj.imaging_data = data;
        obj.neurotar_data = floating;

        obj.importTemplate; % Initial import of the template
    end

    function findPlaceCells(obj)
        obj.makeBins;
        obj.binDFF;
        obj.computeSpatialInformation
        obj.screenPlaceCells;

        plot the stuff
    end

    function makeHeatMaps(obj,cellNum) % i only need this to run on a single cell, which I've specfied?
        obj.makeBins;
        obj.binDFF;
        obj.computeSpatialInformation

        imagesc the hetamaps (I)
end

    function runAnalysis(obj)
        obj.makeBins;
        obj.binDFF;
    end

    function changeTemplate(obj) % This lets you change your template, then re-bin according to the new templaet..
        fprintf('Choose your new template...\n')
        [fn,pn] = uigetfile('.mat');
        obj.importTemplate(pn,fn);
    end

function computeSpatialInformation(obj)
    I = obj.counts / sum(obj.counts(:)) .* obj.binned_DFF    / nanmean(obj.binned_DFF(:)) .* log2(obj.binned_DFF / nanmean(obj.binned_DFF(:)));
I(isnan(I)) = 0; 
obj.spatial_info = sum(I(:)); 
end

function screenPlaceCells(obj)
    if obj.spatial_info > 0
    num_shuffles = 500; 
    shuff_spatial_info = zeros(1, num_shuffles);
    for jj = 1:num_shuffles
        shuff_DFF_binned = zeros(size(DFF_binned)); 
        shuff_start_point = randi(length(cell_responses)); 
        shuff_cell_responses = cell_responses([shuff_start_point:length(cell_responses),...
            1:(shuff_start_point - 1)]); 
        for kk = 1:length(bin_X)
            shuff_DFF_binned(bin_X(kk), bin_Y(kk)) = shuff_DFF_binned(bin_X(kk), bin_Y(kk)) + ...
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

        % Getters and setters
    function out = getBinnedDFF(obj)
            out = obj.DFF_binned;
    end
        
    function out = getTemplate(obj)
            out = obj.template;
    end
        
    function setNumBins(obj,numbins)
            obj.NumBins = numbins;
    end

    function out = exportData(obj) % this packs out all the variables that we need from this object, then we can output things
            out.DFF = obj.binDFF;
            out.counts = obj.counts;
            out.bin_X = obj.bin_X;
            out.bin_Y = obj.bin_Y;
            out.cell_responses = obj.cell_responses;
    end
end
    
methods (Access = protected)
    function importTemplate(obj,pn,fn)
        if nargin < 2
            fprintf('No template provided, loading default empty template... \n')
            obj.template = importdata('D:\Dropbox\Floating Cage\Arena templates\FloatCage_Empty_TEMPLATE.mat');
        else
            obj.template = load([pn fn]);
        end
    end
    
    function makeBins(obj)
        obj.bin_X = obj.template.X(1):(obj.template.X(end) - obj.template.X(1))/(obj.numBins - 1):obj.template.X(end);
        obj.bin_Y = obj.template.Y(1):(obj.template.Y(end) - obj.template.Y(1))/(obj.numBins - 1):obj.template.Y(end);
    end
    
    function binDFF(obj)
        obj.counts = zeros(obj.numBins, obj.numBins);
        dff_temp = zeros(obj.numBins, obj.numBins);
        
        obj.cell_responses = (obj.imaging_data.DFF - min(obj.imaging_data.DFF , [], 2)) ./ ...
            (max(obj.imaging_data.DFF , [], 2) - min(obj.imaging_data.DFF , [], 2)); % normalizing the cell responses to make the heat map more interpretible
        
        for ii = 1:5; %size(cell_responses,1)
            for jj = 1:length(obj.neurotar_data.X)
                x_distance = abs(obj.neurotar_data.X(jj) - bin_X);
                y_distance = abs(obj.neurotar_data.Y(jj) - bin_Y);
                [~, bin_id_X(jj)] = min(x_distance);
                [~, bin_id_Y(jj)] = min(y_distance);
                obj.counts(bin_id_X(jj), bin_id_Y(jj)) = obj.counts(bin_id_X(jj), bin_id_Y(jj)) + 1;
                dff_temp(bin_id_X(jj), bin_id_Y(jj),ii) = dff_temp(bin_id_X(jj), bin_id_Y(jj)) + cell_responses(ii, jj);
            end
        end
        
        obj.DFF_binned = dff_temp ./ obj.counts;
    end

        
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
                        shuff_DFF_binned(bin_X(kk), bin_Y(kk)) = shuff_DFF_binned(bin_X(kk), bin_Y(kk)) + ...
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