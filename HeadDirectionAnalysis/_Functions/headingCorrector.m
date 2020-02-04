function [heading_out] = headingCorrector(decoded_heading, true_heading, bin_width)
	% This is to correct for 180 degree plotting errors, ie this code uses the real heading and chooses the decoded heading closer to that value
	% Written 03Feb2020 KS
	% Updated 

	if nargin < 3 || isempty(bin_width)
		bin_width = 3;
	end
	
	n_bins = 360 / bin_width;
	heading_options = cat(2, decoded_heading, decoded_heading + (n_bins/2));
	heading_options = mod(heading_options, n_bins);

	[~, idx] = min(abs(heading_options - true_heading), [], 2);
	for ii = 1:size(decoded_heading, 1)
		heading_out(ii) = heading_options(ii, idx(ii));
	end

	heading_out = heading_out'; % Transpose to keep the correct dimensions