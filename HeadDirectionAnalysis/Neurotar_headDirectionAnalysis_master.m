% Reformat to MMM###_MM_DD_IMAGING_S.xlsx
clear;

visualize = false;

disp('Choose your matched floating data:')
f_fn = uigetfile('*.mat');
disp('Choose your neural data:')
n_fn = uigetfile('*.mat');

floating = importdata(f_fn);
data = importdata(n_fn);

% Separation happens here
hda = HeadDirectionAnalysis(data, floating);
hda.setHeadingFlag(false)  % Heading false is using alpha
%hda.removeMovingSamples();

%% Working on better analysis
hda.calculateHeadDirectionIdx_direction();
hda.calculatePreferredDirection('vectorsum');  % Max is necessary for the aligned plots

%% Save the data
%hda.saveAnalysisData();
hda.printStats();

%% Beta analyses
vector_length = hda.direction_info(:, 2);


% hda.clusteringAnalysis()
flip_score = hda.calculateFlipScore(); % A measure of how bilobal the things are
histogram(flip_score);
scatter(hda.direction_info(:, 2), flip_score);
ax = gca();
line([0, 0], [ax.YLim(1), ax.YLim(2)])
line([ax.XLim(1), ax.XLim(2)], [0.6, 0.6]) % all numbers from Jacob et al
xlabel('Magnitude')
ylabel('Flip Score')

% 
% % hda.quadrantAnalysisHeadDirection()
% bilobality_info = hda.bilobalityAnalysis(false); % check_flag
% hda.bilobalityPlot();
% 
% sum(hda.is_head_direction & hda.is_bilobed) / sum(hda.is_head_direction)
% 
% save hda_object.mat hda

%% Visualize bilobed cells V2

% Two imagesc plots, showing the arranged HD cells for bilobed and unilobed responses
if visualize
    threshold = 0.5;
    binned_data = hda.getBinnedData();
    bilobed_data = binned_data(flip_score > threshold, :);
    unilobed_data = binned_data(flip_score <= threshold, :);
    
    %sort it
    [~, idx1] =  max(bilobed_data, [], 2);
    [~, sort1] = sort(idx1)
    
    [~, idx2] = max(unilobed_data, [], 2);
    [~, sort2] = sort(idx2);
    
    bilobed_data = bilobed_data(sort1, :);
    unilobed_data = unilobed_data(sort2, :);
    
    subplot(1, 2, 1)
    imagesc(bilobed_data)
    
    subplot(1, 2, 2)
    imagesc(unilobed_data)
end



%% Visualize bilobed cells
if visualize
    binned_data = hda.getBinnedData();
    
    hda.setBinWidth(30);
    wide_plot_data = hda.getBinnedData();
    hda.setBinWidth(3);
    
    is_included = hda.is_head_direction;
    
    for c = 1:size(binned_data, 1)
        if is_included(c)
            subplot(3, 1, 1)
            imagesc(squeeze(plot_data(c, :, :))')
            subplot(3, 1, 2)
            plot(binned_data(c, :))
            subplot(3, 1, 3)
            plot(wide_plot_data(c, :, :))
           % title(bilobality_info.coeffs(c, :))
            pause
        end
    end
end

%% Get the top cells?
if visualize  
    cells2plot = hda.getSignificantCells();
    for c = cells2plot
        subplot(1, 2, 1)
        hda.polarPlot(c);
        subplot(1, 2, 2)
        hda.trajectoryActivityPlot(c);
        pause
        cla
    end
%    hda.alignedPlot();
end


%{







modelfun = @(b,x) b(1) + b(2) * x(:, 1) + ...
    b(3) * exp(-(x(:, 1) - b(4)).^2/b(5)) + b(6) * exp(-(x(:, 1) - b(7)).^2/b(8)); 
[coeffs, r_sqr] = hda.bilobalityFit();

binned_data = hda.getBinnedData();
is_well_fit = r_sqr > 0.5;
well_fit_data = binned_data(r_sqr > 0.5, :);

is_head_direction = hda.getSignificantCells();





x_vals = linspace(1, 12, 100);
for ii = 1:size(binned_data, 1)
    if is_well_fit(ii) & is_head_direction(ii)
    for jj = 1:length(x_vals)
        y(jj) = modelfun(coeffs(ii,:), x_vals(jj));
    end
    plot(x_vals, y)
    hold on
    plot(binned_data(ii, :))
    title([num2str(ii), ' | ' num2str(r_sqr(ii))])
    pause
    hold off
    end
end

%{


[~, mostbilob] = sort(bilobalityIdx);

binnedData = hda.getPlottingData();

prefDir = hda.getPreferredDirection();
[~, sortIdx] = sort(prefDir);
sortedData = binnedData(sortIdx, :);
% 
% rowMin = min(sortedData, [], 2);
% rowMax = max(sortedData, [], 2);
% imagesc(rescale(sortedData, 'InputMin', rowMin, 'InputMax', rowMax));

imagesc(sortedData);


figure
histogram(bilobalityIdx(isHeadDirection))

figure
for ii = 1:10
    subplot(2, 5, ii)
    plot(binnedData(mostbilob(end - (ii - 1)), :))
end



%hda.pointActivityPlot(5);


%{
cells2plot = find(top_cells);
for ii = 1:length(cells2plot)
    subplot(1, 2, 1)
    c = cells2plot(ii);
    plt.polarPlot(c,'LineWidth', 5);
    title(sprintf('Cell #%d', cells2plot(ii)))
    subplot(1, 2, 2)
    imagesc(activity_heatmap(:, :, c))
    axis square
    colorbar
    pause
end


for ii = [18, 176, 404]
    subplot(1, 2, 1)
    plt.polarPlot(ii,'LineWidth', 5);
    title(sprintf('Cell #%d', (ii)))
    subplot(1, 2, 2)
    imagesc(activity_heatmap(:, :, ii))
    axis square
    colorbar
    pause
end
%}


%}

%}