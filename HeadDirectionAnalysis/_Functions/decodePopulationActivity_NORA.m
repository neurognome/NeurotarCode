function [predicted_heading, heading_distribution] = decodePopulationActivity_NORA(tuning_curve, time_series

% Rescale all rows, sets each row to be [0, 1]
tuning_wts = rescale(tuning_curve,...
    'InputMin', min(tuning_curve, [], 2),...
    'InputMax', max(tuning_curve, [], 2));

% cap it
cap = prctile(time_series', 99); % removing mega spikes
for ii = 1:size(time_series, 1)
    temp = time_series(ii, :);
    temp(temp > cap(ii)) = cap(ii);
    time_series(ii, :) = temp;
end

time_series = rescale(time_series,...
    'InputMin', min(time_series, [], 2),...
    'InputMax', max(time_series, [], 2));

pre_post_samples = 10;
heading_distribution = zeros(size(time_series, 2), size(tuning_wts, 2));

for t = 1:size(time_series, 2)
    %temp = (max(mean(timeseries(:, max(t-pre_post_samples, 1) : min(t+pre_post_samples, length(timeseries))), 2), 0) .* tuning_wts) ./ sum(tuning_wts); % nonzero weights are lost
    activity = max(mean(...
        time_series(:, max(t-pre_post_samples, 1) : min(t+pre_post_samples,...
        length(time_series))), 2), 0);
    heading_distribution(t, :) = sum((activity) .* tuning_wts) ./ sum(tuning_wts); % nonzero weights are lost);
end

memory = 5;

for ii = 1:size(heading_distribution, 1)
    current_heading_distribution = mean(heading_distribution(max(1, ii) : ii, :), 1);
    baseline = mean(current_heading_distribution) + 1.3*std(current_heading_distribution);
    
    % break down current vector into logical vector
    curr_peak = current_heading_distribution > baseline;
    
    % eliminate singletons
    for t = 1:size(curr_peak, 2)
        if t == 1
        elseif t == size(curr_peak, 2)
            %t == 2:((length(t))-1)
        else
            if curr_peak(t - 1) == 0 && curr_peak(t + 1) == 0
                curr_peak(t) = 0;
            end
        end
    end
    
    % number of changes in vector/2 to discount negative changes
    change = diff(curr_peak(:));
    num_peak = ceil((nnz(unique(cumsum([true;diff(change)~=0]).*(change~=0))))/2);
    
    %%% Cell array method of peak grouping + meaning
    % find and split peaks
    peak_point = find(curr_peak == 1);
    split = [0 (diff(peak_point))>1];
    split1 = find(split);
    
    
    % group peaks into cell array
    for jj = 1:length(split1) + 1
        
        if num_peak >= 2
            if jj == 1
                peak_group{jj} = peak_point(1:split1(jj) - 1);
            elseif jj == length(split1)+1
                peak_group{jj} = peak_point(split1(jj - 1):end);
            else
                peak_group{jj} = peak_point(split1(jj - 1) : split1(jj) - 1);
            end
            
            avg_cells = cell(1, size(curr_peak, 2));
            % find mean point of each peak
            avg_cells{t} = cellfun(@(x) mean(x(:)), peak_group, 'UniformOutput', false);
            avg_array = vertcat(avg_cells{:});
            avg_peak = cell2mat(avg_array);
        end
    end
    
    % avg peak for <2 peaks
    if num_peak < 2
        avg_peak = mean(peak_point);
    end
    
    if ii <= memory
        predicted_heading(ii) = min(avg_peak);
    else
        previous_heading_distribution = mean(heading_distribution(ii-1, :), 1);
        
        % correlate previous points
        peak_fit = polyfit([(ii-memory):(ii-1)], predicted_heading((ii-memory):(ii-1)), 1);
        if peak_fit(:,1) <= 0.1 && peak_fit(:,1) >= -0.1
            prev_output = ceil(predicted_heading(ii-1));
        else
            prev_output = polyval(peak_fit, ii);
        end
        
        % find current peak closest to the previous peak
        [distance, location] = min(abs(avg_peak - prev_output));
        predicted_heading(ii) = avg_peak(location);
    end
end