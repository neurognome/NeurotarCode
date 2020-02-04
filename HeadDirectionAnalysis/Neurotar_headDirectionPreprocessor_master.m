%Neurotar_headDirectionAnalysis
% Reformat to MMM###_MM_DD_IMAGING_S.xlsx
addpath(genpath('E:\_Code\NeurotarCode'));

clear;

flicker_flag = true;

temp = dir('*.mat');

if isempty(temp)
    error('No matfiles found in current directory, go to the directory containing your matfiles')
end

temp = strcat(temp.name);

if ~contains(temp, 'floating')
    disp('Excel file not extracted yet... doing that now')
    n = NeurotarDataExtractor();
    n.saveData;
end

disp('Choose your floating data:')
f_fn = uigetfile('*.mat');
disp('Choose your neural data:')
n_fn = uigetfile('*.mat');

floating = importdata(f_fn);
data = importdata(n_fn);
hdp = HeadDirectionPreprocessor(data, floating, flicker_flag); % Data, floating, light_dark_flag

hdp.setForceTimeLock(true);
hdp.processData(7); % Input argument is the average window size for meaning across indices
