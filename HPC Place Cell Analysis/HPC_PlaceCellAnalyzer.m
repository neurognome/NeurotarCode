classdef HPC_PlaceCellAnalyzer
properties
end

methods
end
    
end   





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