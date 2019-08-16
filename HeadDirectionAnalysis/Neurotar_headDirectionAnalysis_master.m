%Neurotar_headDirectionAnalysis

% Sample data
%{
load('D:\Code\NeurotarCode\SampleData\floating_data_WTR010_07_24_session1.mat')
load('D:\Code\NeurotarCode\SampleData\TSeries-07242019-1414-001_registered_data.mat')
%}

% if you need it...
%n = NeurotarDataExtractor();

%n.saveData;

hdp = HeadDirectionPreprocessor(data,floating);
%hdp.setForceTimeLock(false);
[data, floating] = hdp.processData;

% Match in python below
hda = HeadDirectionAnalysis(data,floating);

% If you want use raw alpha, but not recommended...
hda.setHeadingFlag(false);

%% Working on better analysis
% Current issue: not sure the best way to be able to figure out if they're actually "head directiony"
[quadrantCorrelations] = hda.calculateHeadDirectionIdx();

isHeadDirection = mean(quadrantCorrelations,2) > 0.2;



%% Visualization
% Visualizing all the cells first...
hda.setData(binned_DFF);


hda.analysisData.exportVar('isDirectionTuned','pref_dir');
cells2plot = find(isHeadDirection);
for ii = 1:length(cells2plot)
    c = cells2plot(ii);
hda.polarPlot(c,'LineWidth',5);
title(sprintf('Cell #%d, prefdir: %0.1f',ii,pref_dir(c)))
pause
end

