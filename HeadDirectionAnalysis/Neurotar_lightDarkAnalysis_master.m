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
    %hda(i_hda).removeMovingSamples(); 
       hda(i_hda).calculatePreferredDirection('vectorsum');
 
    hda(i_hda).calculateHeadDirectionIdx_ori();
    % hda(i_hda).calculateHeadDirectionIdx_direction();
        
  % hda(i_hda).quadrantAnalysisHeadDirection(true);
  % bl_info(i_hda) = hda(i_hda).bilobalityAnalysis(false);
  % flip_score(i_hda, :) = hda(i_hda).calculateFlipScore();
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Above this line is just "preprocessing" of the the light and dark data
% below  this line is the analyses at the current point...

%% You can load previous data for here
lda = LightDarkAnalyzer(hda(1), hda(2));

%lda.correctDarkDrift();

lda.classifyResponses(true); % Check flag

%lda.viewClusteredResponses()

clust_id = lda.getClusters();

%% Looking at each group unaligned
curr_clust = 4
ct = 1;
clear max_idx
for cell_id = find(clust_id == curr_clust)'
    [~, max_idx(ct)] = max(lda.light_data(cell_id ,:));
    ct = ct + 1;
end

group_3 = [lda.light_data(clust_id == curr_clust, :), lda.dark_data(clust_id == curr_clust, :)];
[~, sorting_vec] = sort(max_idx);
group_3 = group_3(sorting_vec, :);
resc = rescale(group_3, 'InputMin', min(group_3, [], 2), 'InputMax', max(group_3, [], 2));
% for ii = 1:size(group_3, 1)
%     resc(ii, :) = zscore(group_3(ii, :));
% end
figure
imagesc(resc)
colormap(flipud(bone))

%% Decoding each group...
tuning_l = hda(1).getBinnedData();
timeseries = [hda(1).getTimeSeries(), hda(2).getTimeSeries()];

for c = unique(clust_id)'

curr_tuning = tuning_l(clust_id == c, :);
curr_ts = timeseries(clust_id == c, :);

decoded_heading = decodePopulationActivity(curr_tuning, curr_ts);

heading = rescale([hda(1).getHeading(); hda(2).getHeading()], 0, 60);


prediction_error = abs(heading - decoded_heading);
prediction_error(prediction_error > 60) = 60;
prediction_error = prediction_error * 360/120;

figure
subplot(3, 2, [1, 2])
plot(heading, 'k:')
hold on
plot(decoded_heading)
hold on
ylabel('heading')


subplot(3, 2, [3, 4])
plot(prediction_error)
title('Prediction Error')
ylabel('Error')

subplot(3, 2, 5)
histogram(prediction_error)
pause
end

%% Gen 2
c = 3
% Decoding each group...
tuning_l = hda(1).getBinnedData();
tuning_d = hda(2).getBinnedData();
timeseries = [hda(1).getTimeSeries() hda(2).getTimeSeries()];

curr_tuning = mean(cat(3, tuning_l(clust_id == c, :), tuning_d(clust_id == c, :)), 3);
curr_ts = timeseries(clust_id == c, :);

[decoded_heading, heading_distribution] = decodePopulationActivity(curr_tuning, curr_ts);

heading = rescale([hda(1).getHeading(); hda(2).getHeading()], 0, 60);


prediction_error = abs(heading - decoded_heading);
prediction_error(prediction_error > 60) = 60;
prediction_error = prediction_error * 360/120;


%% new stufff
rescaled_heading = ...
    rescale(heading_distribution, 'InputMin', min(heading_distribution, [], 2), 'InputMax', max(heading_distribution, [], 2));

clims = [0.5 0.9];
imagesc(movmean(rescaled_heading(10000:12000, :), 1)', clims)
colormap(magma)
axis off


%plot(heading, '--r')
% 
% figure
% subplot(3, 2, [1, 2])
% plot(heading, 'k:')
% hold on
% plot(decoded_heading)
% hold on
% ylabel('heading')
% 
% 
% subplot(3, 2, [3, 4])
% plot(prediction_error)
% title('Predictiol Error')
% ylabel('Error')
% 
% subplot(3, 2, 5)
% histogram(prediction_error)

%%
flip_score = hda(1).calculateFlipScore();

figure
h = histogram(flip_score(clust_id == 1), 10);
hold on
for ii = 2:max(clust_id)
    histogram(flip_score(clust_id == ii), 'BinEdges', h.BinEdges)
end

%{
%% Messing with the decoding now..

% when choosing a tuning curve, should i be using the one from light, dark, or both?

% light first

tuning_l = hda(1).getBinnedData();
timeseries = hda(1).getTimeSeries();

for c = unique(clust_id)'

curr_tuning = tuning_l(clust_id == c, :);
curr_ts = timeseries(clust_id == c, :);

decoded_heading = decodePopulationActivity(curr_tuning, curr_ts);

heading = rescale(hda(1).getHeading(), 0, 120);


prediction_error = abs(heading - predicted_heading);
prediction_error(prediction_error > 60) = 60;
prediction_error = prediction_error * 360/120;

figure
subplot(3, 2, [1, 2])
plot(heading, 'k:')
hold on
plot(predicted_heading)
hold on
ylabel('heading')


subplot(3, 2, [3, 4])
plot(prediction_error)
title('Predictior Error')
ylabel('Error')

subplot(3, 2, 5)
histogram(prediction_error)
pause
end
%}

%

% data_l = hda(1).getBinnedData();
% data_d = hda(2).getBinnedData();
% 
% %% Lobality
% flip_score = hda(1).calculateFlipScore();
% flip_score2 = hda(2).calculateFlipScore();
% flip2 = (flip_score - flip_score2);
% 
% %% tuning
% tuning_light = hda(1).direction_info(:, 2);
% tuning_dark = hda(2).direction_info(:, 2);
% 
% tuning_diff = (tuning_light - tuning_dark);
%
%% Logical skip
%is_stable = true(1, length(hda(1).is_head_direction))';%& hda(2).is_head_direction';
% 
% 
% binned_colors = discretize(flip_score(is_stable), 10);
% color_pool = parula(10);
% 
% colorID = zeros(size(binned_colors, 1), 3);
% for ii = 1:size(binned_colors, 1)
%     colorID(ii, :) = color_pool(binned_colors(ii), :);
% end
% 
% figure
% scatter(tuning_diff(is_stable), flip2(is_stable), [], colorID, 'filled')
% xlabel('tuning difference')
% ylabel('flip diff')
% zlabel('flip score')
% 
% 
% % 
% % hold on
% % y_limits = ylim();
% %line(xlim(), [0.6 0.6], 'Color', [0.7 0.7 0.7], 'LineStyle', '--')
% out = gname;
% 
% temp_l = data_l(is_stable, :);
% temp_d = data_d(is_stable, :);
% 
% cell_idx = str2double(string(cat(1, out.String)));
% 
% light = zeros(length(cell_idx), 120);
% dark = zeros(length(cell_idx), 120);
% ct = 1;
% for c = cell_idx'
%     [~, shift] = max(temp_l(c, :));
%     light(ct, :) = rescale(circshift(temp_l(c, :), 60 - shift));
%     
%     %   [~, shift] = max(data_d(c, :));
% 
%     dark(ct, :) = rescale(circshift(temp_d(c, :), 60 - shift));
%     ct = ct + 1;
% end
% 
% 
% 
% 
% 
% 
% figure
% light = data_l(is_stable, :);
% dark = data_d(is_stable, :);
% for point = out'
%     
%    cell_id = str2double(point.String);
%    
%    plot(light(cell_id, :))
%    hold on
%    plot(dark(cell_id, :));
%    title(cell_id)
%    hold off
%    pause
%    
% end

% figure
% imagesc(out)


% % For checking specific groups later
% rows_to_plot = [210:233];
% figure
% hold on
% plot(curr_data(rows_to_plot, 1:120)', 'Color', [1 0.8 0.8]);
% plot(curr_data(rows_to_plot, 131:end)', 'Color', [0.8 1 0.8]);
% plot(mean(curr_data(rows_to_plot, 1:120)), 'Color', [1, 0, 0], 'LineWidth', 2);
% plot(mean(curr_data(rows_to_plot, 131:end)), 'Color',[0, 1, 0], 'LineWidth', 2);
% 


%data = cat(2, light(is_stable, :), dark(is_stable, :));

% smoothin

% %data = movmean(data, 10, 2);
% methods = {'average',...
%            'centroid',...
%            'complete',...
%            'median',...
%            'single',...
%            'ward',...
%            'weighted'};
% distance = {'euclidean',...
%             'squaredeuclidean',...
%             'seuclidean',...
%             'mahalanobis',...
%             'cityblock',...
%             'minkowski',...
%             'chebychev',...
%             'cosine',...
%             'correlation',...
%             'hamming',...
%             'jaccard',...
%             'spearman'};
% 
% for m = 1:length(methods)
%     for d = 1:length(distance)
% y = pdist(data, distance{d});
% z = linkage(y, methods{m});
% output(m, d) = cophenet(z, y);
%     end
% end

% % My own clustering, round 2 lol
% clusters = unique(t);
% for c = clusters'
%     mean_trace(:, c) = mean(data((t == c), :), 1);
% end
% 
% ct1 = 1;
% for t1 = mean_trace
%     ct2 = 1;
%     for t2 = mean_trace
%         out(ct1, ct2) = corr(t1, t2, 'Type', 'Pearson');
%         ct2 = ct2 + 1;
%     end
%     ct1 = ct1 + 1;
% end
%     
% 
% % time 2 group
% 
% group = zeros(1, size(clusters, 1));
% grp = 1;
% 
% for ii = 1:size(out, 1)
%     to_label = find(out(ii, :) > 0.5);
%     for l = to_label
%     if group(l) == 0
%         group(l) = grp;
%     end
%     end
%     grp = grp  +1;
% end
% 
% % now we cluster the groups...
% tab_grp = tabulate(group);
% sig_groups = tab_grp(tab_grp(:, 2) > 5, 1); % more than 5
% for ii = 1:length(sig_groups)
%     groups_to_find = find(group == sig_groups(ii));
%     running_total = zeros(length(t), 1);
%     for g = groups_to_find
%         running_total = (t == g) +  running_total;
%     end
%     
%     new_groups(logical(running_total)) = ii;
% end

%{
% 
% %% Drift analysis
% light_data = hda(1).getSegmentedData(15); % Bins of 15 leaves it around 1 min per bin
% dark_data = hda(2).getSegmentedData(15);
% 
% fold_flag = false;
% 
% % Calculate the preferred direction in small bins
% for i_cell = 1:size(light_data, 1)
%     for i_segment = 1:size(light_data, 3)
%         light_direction(i_cell, :, i_segment) = hda(1).calculateVectorSum(squeeze(light_data(i_cell, :, i_segment)), fold_flag);
%         dark_direction(i_cell, :, i_segment) = hda(2).calculateVectorSum(squeeze(dark_data(i_cell, :, i_segment)), fold_flag);  
%     end
% end
% 
% light_temp = mean(light_direction, 3);
% dark_temp = mean(dark_direction, 3);
% 
% temp = std(light_direction, [], 3);
% 
% % Stability is the std in the location of the preferred peak over time
% light_stability = 1 ./ temp(:, 1); % inverted so higher values are more stable, just for lookin
% 
% % Tuning is the strength of the tuning
% temp = mean(light_direction, 3);
% light_tuning = temp(:, 2);
% 
% temp = std(dark_direction, [], 3);
% dark_stability = 1 ./ temp(:, 1);
% 
% temp = mean(dark_direction, 3);
% dark_tuning = temp(:, 2);
% 
% tuning_idx = (light_tuning - dark_tuning) ./ (light_tuning + dark_tuning); % focus
% 
% stability_idx = (light_stability - dark_stability) ./ (light_stability + dark_stability);
% 
% 
% % at this point, I think
% is_active = mean(mean(cat(3, data_d, data_l), 3), 2) > 0;

% All the metrics below are for the light conditions
%% Light-Dark Correlation
light_dark_correlation = zeros(1, size(data_l, 1));
for c = 1:size(data_l, 1)
    light_dark_correlation(c) = corr(data_l(c, :)', data_d(c, :)');
end
% Added this scaling factor here in an attempt to incorporate "magnitude" into the correlation
% it's bounded [0, 1]. meaning a very highly dark-responsive neuron 
%scaled_light_dark_correlation = light_dark_correlation .* max(min(mean(data_d, 2) ./ mean(data_l, 2), 1), 0)';
%scaled_light_dark_correlation = scaled_light_dark_correlation';


%% Lobality
flip_score = hda(1).calculateFlipScore();

%% tuning
tuning_light = hda(1).direction_info(:, 2);
tuning_dark = hda(2).direction_info(:, 2);

tuning_diff = tuning_light - tuning_dark;

%% Logical skip
flip_threshold = 0.4;
ldc_threshold = 0.4;
dir_threshold = 0.2;

% %
is_bilobed = flip_score > flip_threshold;
is_lightsensitive = light_dark_correlation < ldc_threshold; % what about negatives here?
%is_direction = directionality > dir_threshold;

% First get rid of cells that aren't even tuned
is_stable = hda(1).is_head_direction';
% is_bilobed = flip_score(is_stable) > flip_threshold;
% is_lightsensitive = scaled_light_dark_correlation(is_stable) < ldc_threshold;
% is_direction = directionality(is_stable) > dir_threshold;

%% Here we can compare stability indices of active bilobed and active unilobed cells
is_heading_cell = ~is_bilobed & is_stable & ~is_lightsensitive';
is_landmark_cell = is_bilobed & is_stable & is_lightsensitive';

% Checks
% fprintf('Number of heading cells: %d; %0.2f%% of the population\n', sum(is_heading_cell), mean(is_heading_cell)*100)
% fprintf('Number of landmark cells: %d; %0.2f%% of the population\n', sum(is_landmark_cell), mean(is_landmark_cell)*100)

% % Comparisons
% % Predictions: 1) Landmark cells should be more stable in light
% 
% fprintf('Heading v Landmark stability: %0.02f vs %0.02f\n',...
% mean(stability_idx(is_heading_cell)),...
% mean(stability_idx(is_landmark_cell)))
% 
% fprintf('Heading v Landmark tuning: %0.02f vs %0.02f\n',...
% mean(tuning_idx(is_heading_cell)),...
% mean(tuning_idx(is_landmark_cell)))


%% Visualize as a sankey heh

% data = [{'all', 'tuned', sum(is_stable)};...
%         {'all', 'untuned', sum(~is_stable)};...
%         {'tuned', 'light sensitive', sum(is_stable & is_lightsensitive')};...
%         {'tuned', 'light insensitive', sum(is_stable & ~is_lightsensitive')};...
%         {'light sensitive', 'visual', sum(is_stable & is_lightsensitive' & is_bilobed)};...
%         {'light sensitive', '???', sum(is_stable & is_lightsensitive' & ~is_bilobed)};...
%         {'light insensitive', 'multimodal', sum(is_stable & ~is_lightsensitive' & is_bilobed)};...
%         {'light insensitive', 'heading', sum(is_stable & ~is_lightsensitive' & ~is_bilobed)}];
% 
% skp = SankeyPlot(data);
% skp.autoRun();

scatter(flip_score(is_stable), light_dark_correlation(is_stable));

figure
for ii = 1:size(data_l, 1)
        %plot(data_l(ii, :))
        if is_stable(ii)
        theta = linspace(0, 2*pi, 121);
        rho = (data_l(ii, :));
        rho = [rho rho(1)];
        polarplot(theta,rho)
        hold on
     %   plot(data_d(ii, :))
             theta = linspace(0, 2*pi, 121);
        rho = (data_d(ii, :));
        rho = [rho rho(1)];
        polarplot(theta, rho)
        title(light_dark_correlation(ii))
        pause
        hold off
        end
end

%{
%% visualize
scatter(hda(1).direction_info(:, 2), flip_score, 'filled')
hold on
% scatter(directionality(~is_lightsensitive), flip_score(~is_lightsensitive), 'filled')

ax = get(gca);
line([dir_threshold, dir_threshold], [ax.YLim(1), ax.YLim(2)])
line([ax.XLim(1), ax.YLim(2)], [flip_threshold, flip_threshold]) % all numbers from Jacob et al


ylabel('Flip score')
xlabel('Directionality')


%{
%% group
is_heading_cell = ~is_bilobed  & is_stable;
is_landmark_cell = is_bilobed  & is_stable;

is_unknown = ~is_bilobed & is_stable;

%% visualize
scatter(directionality(is_lightsensitive), flip_score(is_lightsensitive), 'filled')
hold on
scatter(directionality(~is_lightsensitive), flip_score(~is_lightsensitive), 'filled')

ax = get(gca);
line([dir_threshold, dir_threshold], [ax.YLim(1), ax.YLim(2)])
line([ax.XLim(1), ax.YLim(2)], [flip_threshold, flip_threshold]) % all numbers from Jacob et al


ylabel('Flip score')
xlabel('Directionality')

% This code below is for manually screening cells
for c = 1:size(data_l, 1)
%     if is_multimodal_cell(c)
        y_limit = max([data_l(c, :), data_d(c, :)]);
        subplot(1, 2, 1)
        plot(data_l(c, :))
        ylim([0, y_limit])
        title([num2str(scaled_light_dark_correlation(c)) ' ' num2str(flip_score(c)) ' ' num2str(hda(1).mean_quadrant_correlation(c))])
        % set the ylim so they're on the same scale
        subplot(1, 2, 2)
        plot(data_d(c, :))
        ylim([0, y_limit])
        pause
%     end
end



%{
disp(mean(hda(1).is_head_direction))
disp(mean(hda(2).is_head_direction))

disp(mean(hda(1).is_head_direction & hda(2).is_head_direction))

for c = 1:size(data_l, 1)
    light_dark_correlation(c) = corr(data_l(c, :)', data_d(c, :)');
end
histogram(light_dark_correlation)

% Added this scaling factor here in an attempt to incorporate "magnitude" into the correlation
% it's bounded [0, 1]. meaning a very highly dark-responsive neuron 
scaled_light_dark_correlation = light_dark_correlation .* max(min(mean(data_d, 2) ./ mean(data_l, 2), 1), 0)';

flip_score = hda(1).calculateFlipScore();
directionality = hda(1).direction_info(:, 2);

%{
% Checking amongst the head direction cells, how response are they in each condition?
scatter(mean(data_l(hda(1).is_head_direction, :), 2), mean(data_d(hda(1).is_head_direction, :), 2))
axis square
refline([1, 0])

% Checking the mean quadrant correlation (directional stability)
scatter(hda(1).mean_quadrant_correlation(hda(1).is_head_direction), ...
    hda(2).mean_quadrant_correlation(hda(1).is_head_direction));
refline([1 0])
axis square
%}

%light_dark_idx = hda(2).mean_quadrant_correlation ./ hda(1).mean_quadrant_correlation;
scatter(scaled_light_dark_correlation, flip_score)
%axis([-1 1 -1 1])
ax = gca();
line([0, 0], [ax.YLim(1), ax.YLim(2)])
line([ax.XLim(1), ax.YLim(2)], [0.6, 0.6]) % all numbers from Jacob et al
xlabel('Light Dark Idx')
ylabel('Flip Score')


%% Separating into classes:

% Bearing cells: cells that have consistent response in light or dark
% Landmark cells: cells that are bimodal, but lose their responsiveness in the dark
% Multimodal cells: cells that are bimodal and keep tuning in the dark


%
light_flip = hda(1).calculateFlipScore();
dark_flip = hda(2).calculateFlipScore();

bearing_cells = scaled_light_dark_correlation' > 0.2;
landmark_cells = light_flip > 0.6 & scaled_light_dark_correlation' < 0.2;
multimodal_cells = light_flip > 0.6 & dark_flip > 0.6 & scaled_light_dark_correlation' > 0.2;


scatter(scaled_light_dark_correlation(landmark_cells), light_flip(landmark_cells))
hold on

scatter(scaled_light_dark_correlation(bearing_cells), light_flip(bearing_cells))
scatter(scaled_light_dark_correlation(multimodal_cells), light_flip(multimodal_cells))


% Visualize
for c = 1:size(data_l, 1)
    if hda(1).is_head_direction(c)
        y_limit = max([data_l(c, :), data_d(c, :)]);
        subplot(1, 2, 1)
        plot(data_l(c, :))
        ylim([0, y_limit])
        title(light_dark_correlation(c))
        % set the ylim so they're on the same scale
        subplot(1, 2, 2)
        plot(data_d(c, :))
        ylim([0, y_limit])
        pause
    end
end



%% Drift analysis
for i_hda = 1:2
    hda(i_hda).setBinWidth(3);
end

light_data = hda(1).calculateHeadDirectionIdx(15); % Bins of 15 leaves it around 1 min per bin
dark_data = hda(2).calculateHeadDirectionIdx(15);

[coeffs_l, rsqr_l] = hda.temporaryMegaFit(light_data);
[coeffs_d, rsqr_d] = hda.temporaryMegaFit(dark_data);

is_reliable_l = sum(cat(2, rsqr_l{:}) > 0.6, 2) > 10;
is_reliable_d = sum(cat(2, rsqr_d{:}) > 0.6, 2) > 10;


data = cat(3, light_data, dark_data);

cmap = jet(size(data, 3));
for c = 1:size(data, 1)
    if is_reliable_l(c)
        subplot(3, 1, 1) 
        plot(squeeze(data(c, :, 1)), 'Color', cmap(1, :))
        hold on
        for i_line = 2:size(data, 3)
            plot(squeeze(data(c, :, i_line)), 'Color', cmap(i_line, :))
        end
        hold off
        title(num2str(c))

        subplot(3, 1, 2)
        imagesc(squeeze(data(c, :, :))');
        
        subplot(3, 1, 3)
        plot(mean(data(c, :, :), 3), 'LineWidth', 2)

        pause
    end
end


%{

figure;
for ii = 1:size(light_data, 1)
    if is_reliable_l(ii) || is_reliable_d(ii)
        for jj = 1:4
            subplot(2, 4, jj)
            plot(light_data(ii, :, jj))
            title(['Light: Q' num2str(jj)])
        end
        
        for jj = 1:4
            subplot(2, 4, jj + 4)
            plot(dark_data(ii, :, jj))
            title(['Dark: Q' num2str(jj)])
        end
        pause
    end
end

tabulate(round(bl_info(1).coeffs(hda(1).is_bilobed & hda(1).is_head_direction, 3)))
tabulate(round(bl_info(2).coeffs(hda(2).is_bilobed & hda(2).is_head_direction, 3)))
% Now we compare:

histogram(hda(1).mean_quadrant_correlation)
hold on
histogram(hda(2).mean_quadrant_correlation)

figure; 
plot(hda(1).init_data.frame_F);
hold on
plot(hda(2).init_data.frame_F);
%}

%}

%}
%}
%}