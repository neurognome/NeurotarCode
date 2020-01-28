% Choose a folder containing light dark stuff first
% Go into the folder and run all the analysis pairwise with both, outputting all relevant data nad comparing ya
clear

disp('Choose the folder containing your separated data: ')
base_dir = uigetdir();

cd(base_dir);

files = dir('*.mat');
for i_files = 1:length(files)
    load(files(i_files).name);
end

%% Assign each into object and analyze
hda(1) = HeadDirectionAnalyzer(light_data, light_floating);
hda(2) = HeadDirectionAnalyzer(dark_data, dark_floating);

for i_hda = 1:length(hda)
    hda(i_hda).setHeadingFlag(false);
    hda(i_hda).calculatePreferredDirection('vectorsum');
    hda(i_hda).calculateHeadDirectionIdx_ori();
end

% Getting necessary data
tuning_curve = hda(1).getBinnedData();
timeseries = hda(1).getTimeSeries();


%% Decoding

% Rescale all rows
tuning_wts = rescale(tuning_curve,...
    'InputMin', min(tuning_curve, [], 2),...
    'InputMax', max(tuning_curve, [], 2));

if normalize_flag
    % cap it
    cap = prctile(time_series', 99);
    for ii = 1:size(time_series, 1)
        temp = time_series(ii, :);
        temp(temp > cap(ii)) = cap(ii);
        obj.time_series(ii, :) = temp;
    end
    
    time_series = rescale(time_series,...
        'InputMin', min(time_series, [], 2),...
        'InputMax', max(time_series, [], 2));
end

%             for t = 1:size(timeseries, 2)
%                 temp = (max(timeseries(:, t), 0) .* smoothed_wts) ./ sum(smoothed_wts); % nonzero weights are lost
%                 direction_distributions(t, :) = sum(temp);
%             end

% this is janky af but ok, not sure if i want to keep this
sigmoid_sharpness = 5; % > 5 or else you lose the asymptotes
modelfun = @(x)   (tanh(sigmoid_sharpness * (x - 0.5)) + 1) .* 0.5; %10 defines sharpness

group_method = 'movmean';
switch group_method
    case 'bin'
        bin_size = 10;
        ct=1;
        for t = 1:bin_size:size(time_series, 2)
            temp = (max(mean(time_series(:, t:t + 9), 2), 0) .* tuning_wts) ./ sum(tuning_wts); % nonzero weights are lost
            heading_distribution(ct, :) = sum(temp);
            ct=ct+1;
        end
        
    case 'movmean'
        pre_post_samples = 10;
        heading_distribution = zeros(size(time_series, 2), size(tuning_wts, 2));
        for t = 1:size(time_series, 2)
            %temp = (max(mean(timeseries(:, max(t-pre_post_samples, 1) : min(t+pre_post_samples, length(timeseries))), 2), 0) .* tuning_wts) ./ sum(tuning_wts); % nonzero weights are lost
            activity = max(mean(...
                time_series(:, max(t-pre_post_samples, 1) : min(t+pre_post_samples,...
                length(time_series))), 2), 0);
            heading_distribution(t, :) = sum((activity) .* tuning_wts) ./ sum(tuning_wts); % nonzero weights are lost);
        end
end

% 
% bound = @(x, bl, bu) min(max(x, bl), bu);
% % modelfun = @(b, x) b(1) + ...
% %     max([0, b(2)]) * exp(-(x(:, 1) - bound(b(3), 0, x(end, 1))) .^ 2 / (2 * b(4) .^ 2)) + ...
% %     max([0, b(5)]) * exp(-(x(:, 1) - bound(b(6), 0, x(end, 1))) .^ 2 / (2 * b(4) .^ 2)); % Based on Michael natneuro paper
% 
% modelfun = @(b, x) b(1) + ...
%     max([0, b(2)]) * exp(-(x(:, 1) - bound(b(3), 0, x(end, 1))) .^ 2 / (2 * b(4) .^ 2));
% 
% %% adding some "memory"
% 
% memory = 0;
% for ii = 1:size(heading_distribution, 1)
%     % ok for now, but to be more clever, we should have weighting dependent on the distance to the previous memory, scaling
%     % fctors that decrease with increasing distance...
%     current_heading_distribution = mean(heading_distribution(max(1, ii - memory) : ii, :), 1);
% 
%     coeffs(ii, :) = fitDoubleGaussian(current_heading_distribution, modelfun);
% end
% 
% predicted_heading = coeffs(:, 3); % Location of first (bigger) peak

[~, predicted_heading] = max(heading_distribution, [], 2);

