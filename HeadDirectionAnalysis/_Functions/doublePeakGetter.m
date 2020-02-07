function [heading_out] = doublePeakGetter(heading_distribution, true_heading)

bound = @(x, bl, bu) min(max(x, bl), bu);
% modelfun = @(b, x) b(1) + ...
%     max([0, b(2)]) * exp(-(x(:, 1) - bound(b(3), 0, x(end, 1))) .^ 2 / (2 * b(4) .^ 2)) + ...
%     max([0, b(5)]) * exp(-(x(:, 1) - bound(b(6), 0, x(end, 1))) .^ 2 / (2 * b(4) .^ 2)); % Based on Michael natneuro paper

modelfun = @(b, x) b(1) + ...
    max([0, b(2)]) * exp(-(x(:, 1) - bound(b(3), 0, x(end, 1))) .^ 2 / (2 * b(4) .^ 2));

%% adding some "memory"

memory = 0;
for ii = 1:size(heading_distribution, 1)
    % ok for now, but to be more clever, we should have weighting dependent on the distance to the previous memory, scaling
    % fctors that decrease with increasing distance...
    current_heading_distribution = mean(heading_distribution(max(1, ii - memory) : ii, :), 1);
    coeffs(ii, :) = fitDoubleGaussian(current_heading_distribution, modelfun);
end

predicted_heading = coeffs(:, 3); % Location of first (bigger) peak
end

%% Subfunctions
function [coefficients, fit_quality] = fitDoubleGaussian(data, model)
% The plan is to use a mixture of gaussian models to get a better idea of the bilobality stuff
% Prepare options for linear fitting
opts = statset('nlinfit');
opts.FunValCheck = 'on';
%opts.Display = 'final';

beta0 = prepareBeta0(data);
tbl = table([1:size(data, 2)]', data');
try
    mdl = fitnlm(tbl, model, beta0, 'Options', opts);
    coefficients = mdl.Coefficients{:, 'Estimate'};
    fit_quality = mdl.Rsquared.Adjusted;
catch
    coefficients = NaN;
    fit_quality = 0;
    
    
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
end