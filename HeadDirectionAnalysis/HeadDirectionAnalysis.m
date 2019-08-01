classdef HeadDirectionAnalysis < handle
    
    % Make sure you have the circular statistics toolbox
    properties (Constant = true)
        bin_edges double = [-180:10:180];
        radius    double = [120]; % 120mm
        thresh    double = 0.05; % p-value threshold for significant tuning
    end
    
    properties
        workingData
        analysisData
        plottingData
        
        fixHeadingFlag logical = 1;
    end
    
    properties (Access = protected)
        initData struct
        initFloating struct
    end
    
    methods
        
        function obj = HeadDirectionAnalysis(data,floating)
            obj.initData = data;
            obj.initFloating = floating;
            
            obj.initializeData(data,floating);
        end
        
        function findHeadDirectionCells(obj)
            
            binned_DFF = obj.binDFF(obj.workingData.heading,obj.workingData.DFF);
            obj.analysisData.addData(binned_DFF);
            
            p = obj.detectCells(obj.workingData.DFF);
            obj.analysisData.addData(p);
            
            isDirectionTuned = p < obj.thresh;
            obj.analysisData.addData(isDirectionTuned);       
            
            for c = 1:size(binned_DFF,1)
                if isDirectionTuned(c)
                    obj.polarPlot(c)
                    pause
                end
            end
        end
        
        function setHeadingFlag(obj,val)
            obj.fixHeadingFlag = val;
            obj.initializeData(obj.initData,obj.initFloating);
        end
        
    end
    
    methods (Access = private)
        function initializeData(obj,data,floating)
            DFF = data.DFF;
            X = floating.X;
            Y = floating.Y;
            alpha = floating.alpha;
            phi = floating.phi;
            r = floating.r;
            
            if obj.fixHeadingFlag
                heading = obj.calculateHeading(r,phi,alpha);
            else
                heading = alpha;
            end
            
            obj.workingData  = DataObject('DFF','X','Y','heading','phi','r');
            obj.analysisData = DataObject();
            obj.plottingData = DataObject();
        end
        
        function out = detectCells(obj,data)
            for ii = 1:size(data,1)
                out(ii) = circ_otest(data(ii,:));
            end
        end
        
        function out = binDFF(obj,heading,data)
            groups = discretize(heading,obj.bin_edges);
            u_groups = unique(groups);
            for c = 1:size(data,1)
                for bin = 1:length(u_groups)
                    out(c,bin) = mean(data(c,groups == u_groups(bin)));
                end
            end
        end
        
        function h = calculateHeading(obj,r,phi,alpha)
            h = 2.*asind((r.*sind(alpha - phi))./(2.*obj.radius)) + alpha;
            h = wrapTo180(h); % to account for passing 180 on the thing
        end
        
        function polarPlot(obj,c)
            rho = obj.analysisData.binned_DFF(c,:);
            rho = [rho rho(1)];
            
            theta = obj.bin_edges * pi/180;
            
            polarplot(theta,rho);
        end
        
    end
    
end

