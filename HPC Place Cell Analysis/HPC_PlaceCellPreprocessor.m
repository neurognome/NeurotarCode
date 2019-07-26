classdef HPC_PlaceCellPreprocessor < NeurotarPreProcessor
    

    %% You can add additional proprocessing steps here, but what needs to come out is the raw (relatively) data and the raw (relatively) floating info
    
    properties (Access = protected)
        
    end

    
    methods
        function obj = HPC_PlaceCellPreprocessor(data,floating)
            
            obj@NeurotarPreProcessor(data,floating);           
            
        end
        
        function out = processData(obj)

            processData@NeurotarPreProcessor(obj); % Run the preprocessing steps
    
        end
  
      
end
    
