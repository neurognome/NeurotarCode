function [predicted_heading] = decodePopulationActivity(tuning_curve, time_series, display_flag, normalize_flag)
% At each time point (t), each given neuron contributes a set of vectors, defined by its tuning curve. This set
% of vectors is then scaled by the activity of that neuron at a given time point, and averaged across all neurons
% to get a single population activity vector
if nargin < 3 || isempty(display_flag)
    display_flag = true;
end

if nargin < 4 || isempty(normalize_flag)
    normalize_flag = true;
end

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
    
    obj.time_series = rescale(time_series,...
        'InputMin', min(time_series, [], 2),...
        'InputMax', max(time_series, [], 2));
end

%             for t = 1:size(timeseries, 2)
%                 temp = (max(timeseries(:, t), 0) .* smoothed_wts) ./ sum(smoothed_wts); % nonzero weights are lost
%                 direction_distributions(t, :) = sum(temp);
%             end

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
        bin_size = 1;
        pre_post_samples = 4;
        heading_distribution = zeros(size(time_series, 2), size(tuning_wts, 2));
        for t = 1:size(time_series, 2)
            %temp = (max(mean(timeseries(:, max(t-pre_post_samples, 1) : min(t+pre_post_samples, length(timeseries))), 2), 0) .* tuning_wts) ./ sum(tuning_wts); % nonzero weights are lost
            heading_distribution(t, :) = sum((max(mean(...
                obj.time_series(:, max(t-pre_post_samples, 1) : min(t+pre_post_samples,...
                length(time_series))), 2), 0) .* tuning_wts) ./ sum(tuning_wts)); % nonzero weights are lost);
        end
end

bound = @(x, bl, bu) min(max(x, bl), bu);
modelfun = @(b, x) b(1) + ...
    max([0, b(2)]) * exp(-(x(:, 1) - bound(b(3), 0, x(end, 1))) .^ 2 / (2 * b(4) .^ 2)) + ...
    max([0, b(5)]) * exp(-(x(:, 1) - bound(b(6), 0, x(end, 1))) .^ 2 / (2 * b(4) .^ 2)); % Based on Michael natneuro paper

coeffs = fitDoubleGaussian(heading_distribution, modelfun);

predicted_heading = coeffs(:, 3); % Location of first (bigger) peak


%% Subfunctions
function [coefficients, fit_quality] = fitDoubleGaussian(data, model)
% The plan is to use a mixture of gaussian models to get a better idea of the bilobality stuff
% Prepare options for linear fitting
opts = statset('nlinfit');
opts.FunValCheck = 'on';
%opts.Display = 'final';
coefficients = zeros(size(data, 1), 6);
fit_quality = zeros(size(data, 1), 1);
for c = 1:size(data, 1)
    curr_data = data(c, :);
    
    beta0 = prepareBeta0(curr_data);
    tbl = table([1:size(data, 2)]', data(c, :)');
    try
        mdl = fitnlm(tbl, model, beta0, 'Options', opts);
        coefficients(c, :) = mdl.Coefficients{:, 'Estimate'};
        fit_quality(c) = mdl.Rsquared.Adjusted;
    catch
        coefficients(c, :) = NaN;
        fit_quality(c) = 0;
    end
    
end


function beta0 = prepareBeta0(data)
[sorted_values, sorted_idx] = sort(data, 'descend');
mag1 = max([sorted_values(1), 1]);
loc1 = sorted_idx(1);
std1 = 2;

% find the second peak
wrapN = @(x, N) (1 + mod(x-1, N));
mag2 = mag1 / 2; % Random guess
loc2 = wrapN(loc1 + length(data) / 2, length(data));

beta0 = [0, mag1, loc1, std1, mag2, loc2];


