function [heading_out] = headingCorrector(heading_distribution, true_heading, bin_width)
	% This is to correct for 180 degree plotting errors, ie this code uses the real heading and chooses the decoded heading closer to that value
	% Written 003Feb2020 KS
	% Updated 

	scale = 2; % How many X STD you want
	
	if nargin < 3 || isempty(bin_width)
		bin_width = 3;
	end

	for ii = 1:size(heading_distribution, 1)
		[peak, locs] = findpeaks(heading_distribution(ii, :));
		thresh = mean(peak) + scale * std(peak);
		heading_options = locs(peak > thresh);
		[~, idx] = min(abs(heading_options - true_heading(ii)), [], 2);
		if isempty(heading_options)
			heading_out(ii) = NaN;
		else
			heading_out(ii) = heading_options(idx);
		end
	end

	heading_out = heading_out'; % Transpose to keep the correct dimensions