% Choose a folder containing light dark stuff first
% Go into the folder and run all the analysis pairwise with both, outputting all relevant data nad comparing ya

disp('Choose the folder containing your separated data: ')
base_dir = uigetdir();

cd(base_dir);

files = dir('*.mat');
for i_files = 1:length(files)
    load(files(i_files).name);
end


%% Assign each into object and analyze
hda(1) = HeadDirectionAnalysis(light_data, light_floating);
hda(2) = HeadDirectionAnalysis(dark_data, dark_floating);



for i_hda = 1:length(hda)
    hda(i_hda).setHeadingFlag(false);
    %hda(i_hda).removeMovingSamples(); 
   
    hda(i_hda).calculateHeadDirectionIdx_direction();
    hda(i_hda).calculatePreferredDirection('vectorsum');
        
  % hda(i_hda).quadrantAnalysisHeadDirection(true);
  % bl_info(i_hda) = hda(i_hda).bilobalityAnalysis(false);
  % flip_score(i_hda, :) = hda(i_hda).calculateFlipScore();
end

save light_dark_hda_comparisons.mat hda


%% You can load previous data for here
load('light_dark_hda_comparisons.mat');
data_l = hda(1).getBinnedData();
data_d = hda(2).getBinnedData();

%% Drift analysis
light_data = hda(1).getSegmentedData(15); % Bins of 15 leaves it around 1 min per bin
dark_data = hda(2).getSegmentedData(15);

fold_flag = false;

for i_cell = 1:size(light_data, 1)
    for i_segment = 1:size(light_data, 3)
        light_direction(i_cell, :, i_segment) = hda(1).calculateVectorSum(squeeze(light_data(i_cell, :, i_segment)), fold_flag);
        dark_direction(i_cell, :, i_segment) = hda(2).calculateVectorSum(squeeze(dark_data(i_cell, :, i_segment)), fold_flag);  
    end
end

light_temp = mean(light_direction, 3);
dark_temp = mean(dark_direction, 3);

temp = std(light_direction, [], 3);
light_stability = 1 ./ temp(:, 1); % inverted so higher values are more stable, just for lookin

temp = mean(light_direction, 3);
light_tuning = temp(:, 2);

temp = std(dark_direction, [], 3);
dark_stability = 1 ./ temp(:, 1);

temp = mean(dark_direction, 3);
dark_tuning = temp(:, 2);

tuning_idx = (light_tuning - dark_tuning) ./ (light_tuning + dark_tuning); % focus

stability_idx = (light_stability - dark_stability) ./ (light_stability + dark_stability);


% at this point, I think
is_active = mean(mean(cat(3, data_d, data_l), 3), 2) > 0;

%% All the metrics below are for the light conditions
% %% Light-Dark Correlation
% for c = 1:size(data_l, 1)
%     light_dark_correlation(c) = corr(data_l(c, :)', data_d(c, :)');
% end
% % Added this scaling factor here in an attempt to incorporate "magnitude" into the correlation
% % it's bounded [0, 1]. meaning a very highly dark-responsive neuron 
% scaled_light_dark_correlation = light_dark_correlation .* max(min(mean(data_d, 2) ./ mean(data_l, 2), 1), 0)';
% scaled_light_dark_correlation = scaled_light_dark_correlation';
%% Lobality
flip_score = hda(1).calculateFlipScore();


%% Logical skip
flip_threshold = 0.4;
ldc_threshold = 0.3;
dir_threshold = 0.2;

% %
is_bilobed = flip_score > flip_threshold;
%is_lightsensitive = scaled_light_dark_correlation < ldc_threshold;
%is_direction = directionality > dir_threshold;

% First get rid of cells that aren't even tuned
is_stable = hda(1).is_head_direction';
% is_bilobed = flip_score(is_stable) > flip_threshold;
% is_lightsensitive = scaled_light_dark_correlation(is_stable) < ldc_threshold;
% is_direction = directionality(is_stable) > dir_threshold;

%% Here we can compare stability indices of active bilobed and active unilobed cells
is_heading_cell = ~is_bilobed & is_stable;
is_landmark_cell = is_bilobed & is_stable;

% Checks
fprintf('Number of heading cells: %d; %0.2f%% of the population\n', sum(is_heading_cell), mean(is_heading_cell)*100)
fprintf('Number of landmark cells: %d; %0.2f%% of the population\n', sum(is_landmark_cell), mean(is_landmark_cell)*100)

% Comparisons
% Predictions: 1) Landmark cells should be more stable in light

fprintf('Heading v Landmark stability: %0.02f vs %0.02f\n',...
mean(stability_idx(is_heading_cell)),...
mean(stability_idx(is_landmark_cell)))

fprintf('Heading v Landmark tuning: %0.02f vs %0.02f\n',...
mean(tuning_idx(is_heading_cell)),...
mean(tuning_idx(is_landmark_cell)))


%% visualize
scatter(hda(1).direction_info(:, 2), flip_score, 'filled')
hold on
% scatter(directionality(~is_lightsensitive), flip_score(~is_lightsensitive), 'filled')

ax = get(gca);
line([dir_threshold, dir_threshold], [ax.YLim(1), ax.YLim(2)])
line([ax.XLim(1), ax.YLim(2)], [flip_threshold, flip_threshold]) % all numbers from Jacob et al


ylabel('Flip score')
xlabel('Directionality')

figure

for ii = 1:size(light_data, 1)
        %plot(data_l(ii, :))
        %if is_heading_cell(ii)
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
        title(num2str(ii))
        pause
        hold off
        %end
end


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