function [predicted_heading, heading_distribution] = decodePopulationActivity(tuning_curve, time_series, is_head_direction)

if nargin < 3 || isempty(is_head_direction)
	is_head_direction = true(1, size(tuning_curve, 1));
end
%% Rescale all rows, sets each row to be [0, 1]

%% Updated 25Feb2020 KS
% Changelog: Changed singleton removal to use bwmorph to reduce the chances of getting caps, was getting maxmin'd values a lot
%			 Changed construction of the peak finder for better clarity and performance
%			 Added variable contribution, based on the cell type

%% Parameters
memory = 5; % How many samples to look back for fitting a line for prediction, generally improves performance
pre_post_samples = 10; % how many samples to mean over when creating the heading distribution
% Increasing this value definitely improves performance

tuning_wts = rescale(tuning_curve,...
    'InputMin', min(tuning_curve, [], 2),...
    'InputMax', max(tuning_curve, [], 2));


tuning_wts(~is_head_direction, :) = tuning_wts(~is_head_direction, :) * 0.3; % reduce the contribution of non-hd cells

%{
% cap it
cap = prctile(time_series', 95); % removing mega spikes
for ii = 1:size(time_series, 1)
    temp = time_series(ii, :);
    temp(temp > cap(ii)) = cap(ii);
    time_series(ii, :) = temp;
end


time_series = rescale(time_series,...
    'InputMin', min(time_series, [], 2),...
    'InputMax', max(time_series, [], 2));
%}



%% Get the heading distribution, a matrix of possible headings at each given time point
heading_distribution = zeros(size(time_series, 2), size(tuning_wts, 2));

for t = 1:size(time_series, 2)
    %temp = (max(mean(timeseries(:, max(t-pre_post_samples, 1) : min(t+pre_post_samples, length(timeseries))), 2), 0) .* tuning_wts) ./ sum(tuning_wts); % nonzero weights are lost
    activity = max(mean(...
        time_series(:, max(t-pre_post_samples, 1) : min(t+pre_post_samples,...
        length(time_series))), 2), 0);
    heading_distribution(t, :) = sum((activity) .* tuning_wts) ./ sum(tuning_wts); % nonzero weights are lost);
end


for ii = 1:size(heading_distribution, 1)
    current_heading_distribution = mean(heading_distribution(max(1, ii) : ii, :), 1); % what's the purpose of this line?
    current_heading_distribution = heading_distribution(ii, :); % Does this cover it?
    threshold = mean(current_heading_distribution) + 1 * std(current_heading_distribution);
    
    % break down current vector into logical vector
    curr_peak = current_heading_distribution > threshold;
    
    % eliminate singletons
    % We want to give some leeway to the singletons if they're on the end, but not if they're only on the end, right?

    curr_peak = bwmorph(curr_peak, 'clean'); % this should remove all singletons, but preserve ones on the end, if they have neighbors
    % will have issues if you only have 1 true value on the side, leading to a trail on the other side but probably okay...
   %{
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
%}
    
    % number of changes in vector/2 to discount negative changes
    change = diff(curr_peak(:));
    num_peak = ceil((nnz(unique(cumsum([true;diff(change)~=0]).*(change~=0))))/2); % wtf is going on here... where did you get this
    %% Cell array method of peak grouping + meaning

    % find and split peaks
    peak_point = find(curr_peak == 1);
    split = find([0 (diff(peak_point))>1]);
    
    
    % group peaks into cell array, i want to refactor this later KS 25Feb2020
    if num_peak >= 2
	    for jj = 1:num_peak % A better representation of what we're doing
	    	if num_peak >= 2
	    		if jj == 1
	    			peak_group{jj} = peak_point(1:split(jj) - 1);
	    		elseif jj == length(split)+1
	    			peak_group{jj} = peak_point(split(jj - 1):end);
	    		else
	    			peak_group{jj} = peak_point(split(jj - 1) : split(jj) - 1);
	    		end

	            % avg_cells = cell(1, size(curr_peak, 2));
	            % find mean point of each peak
	            % avg_cells{t} = cellfun(@(x) mean(x(:)), peak_group, 'UniformOutput', false);
	            avg_cells = cellfun(@(x) mean(x(:)), peak_group, 'UniformOutput', false);
	            avg_array = vertcat(avg_cells); % Some of these might not be necessary...
	            avg_peak = cell2mat(avg_array);
	        end
	    end
	    clear peak_group % to ensure that the next run is not affected by previous run
	else
		avg_peak = mean(peak_point);
	end
    
    if ii <= memory
        predicted_heading(ii) = min(avg_peak); % why min?
    else        
        % correlate previous points
        peak_fit = polyfit([(ii-memory):(ii-1)], predicted_heading((ii-memory):(ii-1)), 1);
        if peak_fit(:, 1) <= 0.1 && peak_fit(:, 1) >= -0.1 % If low slope?, ie if slope is [-0.1, 0.1] 
            prev_output = ceil(predicted_heading(ii-1)); % take the highest? why? why not lower?
        else
            prev_output = polyval(peak_fit, ii);
        end
        
        % find current peak closest to the previous peak
        [~, location] = min(abs(avg_peak - prev_output));
        predicted_heading(ii) = avg_peak(location);
    end
end