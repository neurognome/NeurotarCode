%Neurotar_headDirectionAnalysis

% Sample data
% load('D:\Code\NeurotarCode\SampleData\PlaceCells\floating_data_WTR010_07_24_session1.mat')
% load('D:\Code\NeurotarCode\SampleData\PlaceCells\TSeries-07242019-1414-001_registered_data.mat')

load('D:\_HeadDirectionNeurotarData\190731_JBA_MJG022_RSC_150\TSeries-07312019-1437-001\TSeries-07312019-1437-001_registered_data.mat')
load('D:\_HeadDirectionNeurotarData\190731_JBA_MJG022_RSC_150\TSeries-07312019-1437-001\floating_data_MJG022_07_31_session1.mat')

% if you need it...
%n = NeurotarDataExtractor();
%n.saveData;

hdp = HeadDirectionPreprocessor(data,floating);
hdp.setForceTimeLock(true);
[data, floating] = hdp.processData(10); % Input argument is the average window size for meaning across indices


hda = HeadDirectionAnalysis(data, floating);
hda.setHeadingFlag(true)
% If you want use raw alpha, but not recommended...
hda.removeMovingSamples();

%% Working on better analysis
% Current issue: not sure the best way to be able to figure out if they're actually "head directiony"
hda.calculateHeadDirectionIdx();

isHeadDirection = hda.meanQuadrantCorrelation> 0.2;
%disp(num2str(sum(isHeadDirection)));
%disp(num2str(find(isHeadDirection)));

%% Visualization
% Visualizing all the cells first...
binnedData = hda.getPlottingData();
histogram(hda.meanQuadrantCorrelation)


%% Putting prefdirs togehter?
hda.getPreferredDirection('vectorsum')
prefdir = hda.directionInfo(:, 1);

for c = 1:size(binnedData)
    rotated(c, :) = circshift(binnedData(c, :), 9 - prefdir(c));
end

plot(mean(rotated(isHeadDirection, :)))

%% Get the top cells?
top_cells = hda.meanQuadrantCorrelation > prctile(hda.meanQuadrantCorrelation, 95);

plt = RawDataPlots();

binnedData = hda.getPlottingData();

plt.setData(binnedData);

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
