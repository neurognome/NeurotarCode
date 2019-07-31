classdef PlaceCellAnalysis < NeurotarPreProcessor
    
    properties (Constant = true)
        n_bins = 10;
    end
    
    properties (SetAccess = protected)
        Template
        bin_X
        bin_Y
        cell_responses
        DFF_binned
        counts
    end
    
    methods
        function obj = PlaceCellAnalysis(data,floating)
            obj@NeurotarPreProcessor(data,floating)
            
            obj.Template = importdata('D:\Dropbox\Floating Cage\Arena templates\FloatCage_Empty_TEMPLATE.mat');
            
            obj.bin_X = obj.Template.X(1):(obj.Template.X(end) - obj.Template.X(1))/(obj.n_bins - 1):obj.Template.X(end);
            obj.bin_Y = obj.Template.Y(1):(obj.Template.Y(end) - obj.Template.Y(1))/(obj.n_bins - 1):obj.Template.Y(end);
            
            obj.cell_responses = (obj.imaging_data.DFF - min(obj.imaging_data.DFF , [], 2)) ./ ...
                (max(obj.imaging_data.DFF , [], 2) - min(obj.imaging_data.DFF , [], 2)); % normalizing the cell responses to make the heat map more interpretible
            
            for c = 1:1;%size(obj.imaging_data.DFF,1)
            [bin_idx_X,bin_idx_Y] = obj.hpcSpatialPlotter(c);
            
            [I,place_cell,spatial_info] = obj.definingPlaceCells(obj.cell_responses(c,:),bin_idx_X,bin_idx_Y);
            end
        end  
    end
    
    
    
    
    methods (Access = private)
        
        function [bin_id_X,bin_id_Y] = hpcSpatialPlotter(obj,ct)
                obj.counts = zeros(obj.n_bins, obj.n_bins); 
                obj.DFF_binned = zeros(obj.n_bins, obj.n_bins); 

            for jj = 1:length(obj.neurotar_data.X)
                x_distance = abs(obj.neurotar_data.X(jj) - obj.bin_X);
                y_distance = abs(obj.neurotar_data.Y(jj) - obj.bin_Y);
                [~, bin_id_X(jj)] = min(x_distance);
                [~, bin_id_Y(jj)] = min(y_distance);
                obj.counts(bin_id_X(jj), bin_id_Y(jj)) = obj.counts(bin_id_X(jj), bin_id_Y(jj)) + 1;
                obj.DFF_binned(bin_id_X(jj), bin_id_Y(jj)) = obj.DFF_binned(bin_id_X(jj), bin_id_Y(jj)) + obj.cell_responses(ct, jj);
            end
            
            
            obj.DFF_binned = obj.DFF_binned ./ obj.counts;
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

%{
function parietalFiringRatePlotter(obj)
            N = size(obj.imaging_data.DFF, 1);
            for ii = 1:N
                % this loop is a bit redundant (since it finds the counts each time)
                % but it doesn't cost much time and it makes things tidier
                counts = zeros(obj.n_bins, obj.n_bins);
                DFF_binned = zeros(obj.n_bins, obj.n_bins);
                distance2objects_bins = 0:(max([max(obj.bin_X) - min(obj.bin_X), max(obj.bin_Y) - min(obj.bin_Y)]) / 49):...
                    max([max(obj.bin_X) - min(obj.bin_X), max(obj.bin_Y) - min(obj.bin_Y)]);
                orientation_bins = -100:10:100;
                DFF_circles = zeros(length(obj.Template.circles), length(distance2objects_bins), length(orientation_bins));
                circle_counts = zeros(length(obj.Template.circles), length(distance2objects_bins), length(orientation_bins));
                DFF_polygons = zeros(length(obj.Template.polygons), length(distance2objects_bins), length(orientation_bins));
                polygon_counts = zeros(length(obj.Template.polygons), length(distance2objects_bins), length(orientation_bins));
                circle_centroids = zeros(length(obj.Template.circles), 2);
                polygon_centroids = zeros(length(obj.Template.polygons), 2);
                
                for jj = 1:length(obj.Template.circles)
                    cent = regionprops(obj.Template.circles{jj}, 'centroid');
                    cent = cell2mat(struct2cell(cent));
                    circle_centroids(jj, :) = [obj.bin_X(round(cent(1) * n_bins / 100)), ...
                        obj.bin_Y(round(cent(2) * obj.n_bins / 100))];
                end
                for jj = 1:length(obj.Template.polygons)
                    cent = regionprops(obj.Template.polygons{jj}, 'centroid');
                    cent = cell2mat(struct2cell(cent));
                    polygon_centroids(jj, :) = [obj.bin_X(round(cent(1) * obj.n_bins / 100)), ...
                        obj.bin_Y(round(cent(2) * obj.n_bins / 100))];
                end
                
                bin_id_X = zeros(1, length(obj.neurotar_data.X));
                bin_id_Y = zeros(1, length(obj.neurotar_data.Y));
                for jj = 1:length(obj.neurotar_data.X)
                    x_distance = abs(obj.neurotar_data.X(jj) - obj.bin_X);
                    y_distance = abs(obj.neurotar_data.Y(jj) - obj.bin_Y);
                    [~, bin_id_X(jj)] = min(x_distance);
                    [~, bin_id_Y(jj)] = min(y_distance);
                    counts(bin_id_X(jj), bin_id_Y(jj)) = counts(bin_id_X(jj), bin_id_Y(jj)) + 1;
                    DFF_binned(bin_id_X(jj), bin_id_Y(jj)) = DFF_binned(bin_id_X(jj), bin_id_Y(jj)) + obj.cell_responses(ii, jj);
                    
                    for kk = 1:length(obj.Template.circles)
                        distance = sqrt((obj.neurotar_data.X(jj) - circle_centroids(kk, 1)) ^ 2 ...
                            + (obj.neurotar_data.Y(jj) - circle_centroids(kk, 2)) ^ 2);
                        [~, bin_id_distance] = min(abs(distance - distance2objects_bins));
                        theta = 180 * asin((obj.neurotar_data.X(jj) - circle_centroids(kk, 1)) / distance) / pi;
                        [~, bin_id_orientation] = min(abs(theta - orientation_bins));
                        DFF_circles(kk, bin_id_distance, bin_id_orientation) = ...
                            DFF_circles(kk, bin_id_distance, bin_id_orientation) + ...
                            obj.cell_responses(ii, jj);
                        circle_counts(kk, bin_id_distance, bin_id_orientation) = ...
                            circle_counts(kk, bin_id_distance, bin_id_orientation) + 1;
                    end
                    
                    for kk = 1:length(obj.Template.polygons)
                        distance = sqrt((obj.neurotar_data.X(jj) - polygon_centroids(kk, 1)) ^ 2 ...
                            + (obj.neurotar_data.Y(jj) - polygon_centroids(kk, 2)) ^ 2);
                        [~, bin_id_distance] = min(abs(distance - distance2objects_bins));
                        theta = 180 * asin((obj.neurotar_data.X(jj) - circle_centroids(kk, 1)) / distance) / pi;
                        [~, bin_id_orientation] = min(abs(theta - orientation_bins));
                        DFF_polygons(kk, bin_id_distance, bin_id_orientation) = ...
                            DFF_polygons(kk, bin_id_distance, bin_id_orientation) + ...
                            obj.cell_responses(ii, jj);
                        polygon_counts(kk, bin_id_distance, bin_id_orientation) = ...
                            polygon_counts(kk, bin_id_distance, bin_id_orientation) + 1;
                    end
                end
            end
        end
    
%}
    
    
%}
