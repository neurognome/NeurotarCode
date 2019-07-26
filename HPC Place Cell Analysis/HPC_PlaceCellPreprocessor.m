classdef HPC_PlaceCellPreprocessor < NeurotarPreProcessor
    
    
    
    properties (Access = protected)
        template struct
        
        DFF_binned double      
        counts double
        cell_responses double
        bin_x double
        bin_y double
    
    
        NumBins double = 10;
    end
    
    methods
        function obj = HPC_PlaceCellPreprocessor(data,floating)
            
            obj@NeurotarPreProcessor(data,floating);           
            
            obj.importTemplate; % Initial import of the template
        end
        
        function out = processData(obj)
            processData@NeurotarPreProcessor(obj); % Run the preprocessing steps
            
            [bin_X,bin_Y]=obj.makeBins;
            [counts,DFF_binned,bin_X,bin_Y,cell_responses] = obj.binDFF(obj,bin_X,bin_Y);
        end
        
        function out = updateData(obj) % Updating, doesn't include initial processing
            [bin_X,bin_Y] = obj.makeBins;
            [counts,DFF_binned,bin_X,bin_Y,cell_responses] = obj.binDFF(obj,bin_X,bin_Y);
            
           %% last step, find the cells
        end
        
        function changeTemplate(obj) % This lets you change your template, then re-bin according to the new templaet..
            fprintf('Choose your new template...\n')
            [fn,pn] = uigetfile('.mat');
            obj.importTemplate(pn,fn);
        end
        
        function out = getBinnedDFF(obj)
            out = obj.DFF_binned;
        end
        
        function out = getTemplate(obj)
            out = obj.template;
        end
        
        function setNumBins(obj,numbins)
            obj.NumBins = numbins;
        end
            
    end
    
    methods (Access = public)
        function importTemplate(obj,pn,fn)
            if nargin < 2
                fprintf('No template provided, loading default empty template... \n')
                obj.template = importdata('D:\Dropbox\Floating Cage\Arena templates\FloatCage_Empty_TEMPLATE.mat');
            else
                obj.template = load([pn fn]);
            end
        end
        
        function [out1,out2] = makeBins(obj)
            
            bin_X = obj.template.X(1):(obj.template.X(end) - obj.template.X(1))/(obj.NumBins - 1):obj.template.X(end);
            bin_Y = obj.template.Y(1):(obj.template.Y(end) - obj.template.Y(1))/(obj.NumBins - 1):obj.template.Y(end);
            out1 = bin_X;
            out2 = bin_Y;
        end
        
        function [counts,DFF_binned,bin_X,bin_Y,cell_responses] = binDFF(obj,bin_X,bin_Y)
            counts = zeros(obj.NumBins, obj.NumBins);
            DFF_binned = zeros(obj.NumBins, obj.NumBins);
            
            cell_responses = (obj.imaging_data.DFF - min(obj.imaging_data.DFF , [], 2)) ./ ...
                (max(obj.imaging_data.DFF , [], 2) - min(obj.imaging_data.DFF , [], 2)); % normalizing the cell responses to make the heat map more interpretible
            
            for ii = 1:5; %size(cell_responses,1)
                for jj = 1:length(obj.neurotar_data.X)
                    x_distance = abs(obj.neurotar_data.X(jj) - bin_X);
                    y_distance = abs(obj.neurotar_data.Y(jj) - bin_Y);
                    [~, bin_id_X(jj)] = min(x_distance);
                    [~, bin_id_Y(jj)] = min(y_distance);
                    counts(bin_id_X(jj), bin_id_Y(jj)) = counts(bin_id_X(jj), bin_id_Y(jj)) + 1;
                    DFF_binned(bin_id_X(jj), bin_id_Y(jj),ii) = DFF_binned(bin_id_X(jj), bin_id_Y(jj)) + cell_responses(ii, jj);
                end
            end
            
            DFF_binned = DFF_binned ./ counts;
        end

            
    end
end
    
