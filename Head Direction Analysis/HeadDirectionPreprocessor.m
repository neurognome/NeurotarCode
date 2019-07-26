classdef HeadDirectionPreprocessor < NeurotarPreProcessor
    properties (Constant = true)
        bin_edges = [-180:10:180];
    end
    
    properties (Access = protected)
        processed_DFF
    end
    
    
    methods
        function obj = HeadDirectionPreprocessor(data,floating,varargin)
            obj@NeurotarPreProcessor(data,floating,varargin{:});
            
            dff = obj.getImagingData.DFF;
            direction = obj.getNeurotarData.alpha;
            
            obj.processed_DFF = zeros(size(dff,1),length(obj.bin_edges)-1); % initialize the processedDFF matrix
            groups = discretize(direction,obj.bin_edges);                   % discretize the directions
            
            obj.binResponse(dff,groups)                                     % Initial binning, always done
        end
        
        function smoothResponse(obj,smooth_factor)
            if nargin == 0
                smooth_factor = round(length(obj.processed_DFF)/10);
            end
            obj.processed_DFF = smoothdata(obj.processed_DFF,2,'movmean',smooth_factor);
        end
        
        % Getters
        function data = getProcessedDFF(obj)
            data = obj.processed_DFF;
        end
    end
    
    
    
    methods (Access = private)
        function binResponse(obj,dff,groups)
            u_groups = unique(groups);
            for c = 1:size(dff,1)
                for bin = 1:length(u_groups)
                    obj.processed_DFF(c,bin) = mean(dff(c,groups == u_groups(bin)));
                end
            end
        end
    end
    
end

