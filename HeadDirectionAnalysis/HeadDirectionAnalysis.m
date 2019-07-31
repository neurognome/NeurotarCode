classdef HeadDirectionAnalysis < handle
       properties (Constant = true)
        bin_edges double = [-180:10:180];
        radius    double = [60]; % 60mm?
    end
    
    properties 
        workingData
        analysisData
        plottingData
        
        initData struct
        initFloating struct
    end 
    
    methods
        function obj = HeadDirectionAnalysis(data,floating)
            obj.initData = data;
            obj.initFloating = floating;
            
            DFF = data.DFF;
            X = floating.X;
            Y = floating.Y;
            alpha = floating.alpha;
            obj.workingData  = dataObject(DFF,X,Y,alpha);
            obj.analysisData = dataObject();
            obj.plottingData = dataObject();
        end
        
        function findHeadDirectionCells();
            calculateHeading
            binResponse
        end
        
        function [x,y] = calculateHeading(obj)
            heading = 2*asind(()
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

