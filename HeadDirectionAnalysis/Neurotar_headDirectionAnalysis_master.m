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
%hdp.setForceTimeLock(true);
%hdp.processData;
hda = HeadDirectionAnalysis(hdp.workingData.get('data'), hdp.workingData.get('floating'));



% If you want use raw alpha, but not recommended...
hda.setHeadingFlag(false);


%% Preprocessing 
binned_DFF = hda.binDFF(hda.workingData.get('heading'),hda.workingData.get('DFF')); % bin DFF
DFF = hda.workingData.get('DFF');
% Not currently running below because I don't have the toolbox
%isDirectionTuned = hda.detectCells(DFF); % Detect non-circularly uniform cells
pref_dir = hda.getPreferredDirection(binned_DFF); % Get preferred directions

hda.analysisData.add('binned_DFF','pref_dir');



%% Working on better analysis
% Current issue: not sure the best way to be able to figure out if they're actually "head directiony"
[headDirIdx] = hda.calculateHeadDirectionIdx();

isHeadDirection = mean(headDirIdx,2) > 0.2;

% what's the best way to determine that? Idk we'll work on this...


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



