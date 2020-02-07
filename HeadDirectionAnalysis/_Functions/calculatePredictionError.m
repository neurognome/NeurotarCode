function prediction_error = calculatePredictionError(decoded_heading, true_heading)
%{
	prediction_error = abs(true_heading - decoded_heading);
    prediction_error(prediction_error >  ) = 60;
    prediction_error = prediction_error * 360/120;
%}

	for ii = 1:length(decoded_heading)
		if decoded_heading(ii) < true_heading(ii)
			t = true_heading(ii);
			true_heading(ii) = decoded_heading(ii);
			decoded_heading(ii) = t;
		end

		dist = decoded_heading(ii) - true_heading(ii) + 1;
		dist2 = true_heading(ii) - 1 + 360 - decoded_heading(ii);
		prediction_error(ii) = min(dist, dist2);
	end
end