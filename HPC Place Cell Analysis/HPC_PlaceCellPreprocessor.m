classdef HPC_PlaceCellPreprocessor < NeurotarPreProcessor
    
    
    properties (Constant = true)
        n_bins = 10;
    end
    
    properties (Access = public)
        template
        DFF_binned
    end
    
    methods
        function obj = HPC_PlaceCellPreprocessor(data,floating,varargin)
            obj@NeurotarPreProcessor(data,floating,varargin{:});
            
            obj.importTemplate;
            
            obj.DFF_binned = obj.dffBinner;
        end
        
        function changeTemplate(obj) % This lets you change your template, then re-bin according to the new templaet..
            fprintf('Choose your new template...\n')
            [fn,pn] = uigetfile('.mat');
            obj.importTemplate(pn,fn);
            fprintf('Re-binning DFF..');
            obj.DFF_binned = obj.dffBinner;
        end
        
        function out = getBinnedDFF
            
        end
        
    end
    methods (Access = private)
        function importTemplate(obj,pn,fn)
            if nargin < 2
                fprintf('No template provided, loading default empty template... \n')
                obj.template = importdata('D:\Dropbox\Floating Cage\Arena templates\FloatCage_Empty_TEMPLATE.mat');
            else
                obj.template = load([pn fn]);
            end
        end
        
        
        function out = dffBinner(obj)
            counts = zeros(obj.n_bins, obj.n_bins);
            dff_temp = zeros(obj.n_bins, obj.n_bins);
            
            bin_X = obj.template.X(1):(obj.template.X(end) - obj.template.X(1))/(obj.n_bins - 1):obj.template.X(end);
            bin_Y = obj.template.Y(1):(obj.template.Y(end) - obj.template.Y(1))/(obj.n_bins - 1):obj.template.Y(end);
            
            cell_responses = (obj.imaging_data.DFF - min(obj.imaging_data.DFF , [], 2)) ./ ...
                (max(obj.imaging_data.DFF , [], 2) - min(obj.imaging_data.DFF , [], 2)); % normalizing the cell responses to make the heat map more interpretible
            for ii = 1:5; %size(cell_responses,1)
                for jj = 1:length(obj.neurotar_data.X)
                    x_distance = abs(obj.neurotar_data.X(jj) - bin_X);
                    y_distance = abs(obj.neurotar_data.Y(jj) - bin_Y);
                    [~, bin_id_X(jj)] = min(x_distance);
                    [~, bin_id_Y(jj)] = min(y_distance);
                    counts(bin_id_X(jj), bin_id_Y(jj)) = counts(bin_id_X(jj), bin_id_Y(jj)) + 1;
                    dff_temp(bin_id_X(jj), bin_id_Y(jj),ii) = dff_temp(bin_id_X(jj), bin_id_Y(jj)) + cell_responses(ii, jj);
                end
            end
            
            out = dff_temp ./ counts;
        end
        
    end
    
end
end
