classdef HeadDirectionAnalysis < RawDataPlots & SummaryDataPlots
    
    % Make sure you have the circular statistics toolbox
    properties (Constant = true)
        bin_width double = 10;
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

        function findHeadDirectionCells(obj) % just a wrapper for full analysis later
        end
%% Setters and Getters       
        function setHeadingFlag(obj,val)
            obj.fixHeadingFlag = val;
            obj.initializeData(obj.initData,obj.initFloating);
        end
    end
    
    
    methods (Access = public)
        function initializeData(obj,data,floating)
            try % temporarily to account of spikeInferenc/noSpikeInference
                DFF = data.spikes;
            catch
                fprintf('Not spike inferred, using DFF\n')
                DFF = data.DFF;
            end
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
                otest(ii) = circ_otest(data(ii,:)); %nonuniformity
                if otest(ii) < obj.thresh
                    rtest(ii) = circ_rtest(data(ii,:)); %unimodality
                else 
                    rtest(ii) = 1;
                end
            end
            out = (otest<obj.thresh) & (rtest<obj.thresh);
        end
        
        function out = getPreferredDirection(obj,data)
            fprintf('Getting preferred directions \n');
           for ii = 1:size(data,1)
              % lm = fit([1:length(data(ii,:))]',data(ii,:)','gauss1','Upper',[Inf 36 Inf], 'Lower', [-Inf 0 -Inf]);
             %  out(ii) = lm.b1;
             [~,out(ii)] = max(data(ii,:)); 
           end 
           
           out=out*obj.bin_width - obj.bin_width; % convert to degrees, accounting for the extra 10 you get because -180:10:180
        end
        function out = binDFF(obj,heading,data)
            bin_edges = -180:obj.bin_width:180; % Because alpha is [-180,180];
            groups = discretize(heading,bin_edges);
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

    end
    
end

        
        %% Major analysis code 
        
        % when you're done messing around in the script, you can incorporate as a big thing, but until then, run the individual functions in the script ok
        
        %{
        function findHeadDirectionCells(obj)
            
            binned_DFF = obj.binDFF(obj.workingData.get('heading'),obj.workingData.get('DFF'));
            obj.analysisData.add(binned_DFF);
            
            p = obj.detectCells(obj.workingData.get('DFF'));
            pref_dir = obj.getPreferredDirection(obj.analysisData.get('binned_DFF'));
            
            obj.analysisData.add('p','pref_dir');
            
            isDirectionTuned = p < obj.thresh;
            obj.analysisData.add(isDirectionTuned);

        end
        %}

