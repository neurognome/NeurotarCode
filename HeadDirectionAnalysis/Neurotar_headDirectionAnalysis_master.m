%Neurotar_headDirectionAnalysis

load('D:\Code\NeurotarCode\SampleData\floating_data_WTR010_07_24_session1.mat')
load('D:\Code\NeurotarCode\SampleData\TSeries-07242019-1414-001_registered_data.mat')

hdp = HeadDirectionPreprocessor(data,floating);
hdp.setForceTimeLock(true);
hdp.processData;

hda = HeadDirectionAnalysis(hdp.workingData.data,hdp.workingData.floating);

% If you want use raw alpha, but not recommended...
%hda.setHeadingFlag(false);

hda.findHeadDirectionCells();

hda.analysisData.exportToVar('binned_DFF');


circ_mean(binned_DFF')

