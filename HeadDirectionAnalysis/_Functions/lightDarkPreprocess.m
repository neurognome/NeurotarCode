function hda = lightDarkPreprocess()

addpath(genpath('E:\_Code\NeurotarCode\HeadDirectionAnalysis'))

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
    hda(i_hda).calculateHeadDirection();
end