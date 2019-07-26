classdef HeadDirectionVisualize < handle
    
    properties (Access = protected)
        data 
    end
    
    methods
        function obj = HeadDirectionVisualize(data)
           obj.data = data; 
        end
        
        
          function polarPlot(obj,cellNum,rectify_flag)
            rho = obj.data(cellNum,:);
            
            theta = [0 :(2*pi)/length(rho): 2*pi];
            
            rho = [rho rho(1)];
            if rectify_flag
                rho = rho - min(rho(:));
            end
            polarplot(theta,rho)
        end
        
    end
end
