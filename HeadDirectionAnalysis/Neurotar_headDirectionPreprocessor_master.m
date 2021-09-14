%Neurotar_headDirectionAnalysis
% Reformat to MMM###_MM_DD_IMAGING_S.xlsx
addpath(genpath('E:\_Code\NeurotarCode'));
addpath(genpath('E:\_Code\Other Helper Functions (not from me)\TDMS_Reader'));
% clear;

flicker_flag = false;

temp = dir('*_registered_data.mat');

if isempty(temp)
    error('No matfiles found in current directory, go to the directory containing your matfiles')
end

temp = strcat(temp.name);

if ~contains(temp, '_raw_stimulus.mat')  
    disp('TDMS file not extracted yet... doing that now')
    floating = convertTDMS(false);
end

% disp('Choose your neural data:')
% n_fn = uigetfile('*.mat');
n_fn = temp;

data = importdata(n_fn);
npp = NeurotarPreProcessor(data, floating); % Data, floating, light_dark_flag

npp.setForceTimeLock(true);
npp.processData(7); % Input argument is the average window size for meaning across indices
% delete(n_fn);
