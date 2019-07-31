classdef HeadDirectionPreprocessor < NeurotarPreProcessor
    properties (Access = protected)
    end
    
    
    methods
        function obj = HeadDirectionPreprocessor(data,floating,varargin)
            obj@NeurotarPreProcessor(data,floating,varargin{:});
        end 
    end
    
    
    
end

