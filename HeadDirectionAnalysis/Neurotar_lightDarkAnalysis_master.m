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
    hda(i_hda).removeMovingSamples();
    
    hda(i_hda).calculateHeadDirectionIdx();
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


%% All the metrics below are for the light conditions
%% HD-ness
directionality = rescale(hda(1).direction_info(:, 2));

%histogram(hda(1).direction_info(:, 1))

%% Light-Dark Correlation
for c = 1:size(data_l, 1)
    light_dark_correlation(c) = corr(data_l(c, :)', data_d(c, :)');
end
% Added this scaling factor here in an attempt to incorporate "magnitude" into the correlation
% it's bounded [0, 1]. meaning a very highly dark-responsive neuron 
scaled_light_dark_correlation = light_dark_correlation .* max(min(mean(data_d, 2) ./ mean(data_l, 2), 1), 0)';
scaled_light_dark_correlation = scaled_light_dark_correlation';
%% Lobality
flip_score = hda(1).calculateFlipScore();


%% Logical skip
flip_threshold = 0.6;
ldc_threshold = 0.3;
dir_threshold = 0.2;

%
is_bilobed = flip_score > flip_threshold;
is_lightsensitive = scaled_light_dark_correlation < ldc_threshold;
is_direction = directionality > dir_threshold;

mean(is_lightsensitive & (is_bilobed & is_direction)) % How many cells are bilobed and direction tuned
is_lightsensitive & (~is_bilobed & is_direction) % How many cells are bilobed and direction tuned


%% visualize
scatter(directionality(is_lightsensitive), flip_score(is_lightsensitive), 'filled')
hold on
scatter(directionality(~is_lightsensitive), flip_score(~is_lightsensitive), 'filled')

ax = get(gca);
line([dir_threshold, dir_threshold], [ax.YLim(1), ax.YLim(2)])
line([ax.XLim(1), ax.YLim(2)], [flip_threshold, flip_threshold]) % all numbers from Jacob et al


ylabel('Flip score')
xlabel('Directionality')

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