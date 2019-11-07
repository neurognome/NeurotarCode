function [light_data, dark_data, light_floating, dark_floating] = lightDarkSeparator(first_segment, floating, data, save_flag)
% This function is used to separate the light and dark responses into two separate data and floating files. These can then
% individually be run through the head direction analysis pipeline to get differences between light and darkness. you shuold
% run this after the full 2p processing pipeline

%% Parameters
segment_duration = 300; % How long each segment is in seconds
fs = 10;

%% Defaults
if nargin < 2
    fprintf('Choose your floating:\n')
    load(uigetfile('.mat'));
end

if nargin < 3
    fprintf('Choose your data:\n')
    load(uigetfile('.mat'))
end


if nargin < 4
    save_flag = true;
end
segment_length = segment_duration * fs;
recording_length = size(data.DFF, 2);

switch first_segment
    case 'light'
        [light_data, dark_data] = separateLightAndDark(data, segment_length, recording_length);
        [light_floating, dark_floating] = separateLightAndDark(floating, segment_length, recording_length);
    case 'dark'
        [dark_data, light_data] = separateLightAndDark(data, segment_length, recording_length);
        [dark_floating, light_floating] = separateLightAndDark(floating, segment_length, recording_length);
end

if save_flag
    mkdir('Separated Data')
    old = cd('Separated Data');
    
    save light_data.mat light_data
    save dark_data.mat dark_data
    save light_floating.mat light_floating
    save dark_floating.mat dark_floating
    
    cd(old);
end
end

%% Subfunctions

function [data1, data2] = separateLightAndDark(data, segment_length, recording_length)
%% Get all the fields in the data structure

data1 = struct();
data2 = struct();
struct_fields = fields(data);
for i_field = 1:length(struct_fields)
    if ismember(recording_length, size(data.(struct_fields{i_field})))
        [data1.(struct_fields{i_field}), data2.(struct_fields{i_field})] = ...
            separateData(data.(struct_fields{i_field}), segment_length, recording_length);
    else
        data1.(struct_fields{i_field}) = data.(struct_fields{i_field});
        data2.(struct_fields{i_field}) = data.(struct_fields{i_field});
    end
end
end

function [data1, data2] = separateData(data, segment_length, recording_length)
% For separating out a single data

n_segments = recording_length / segment_length;

if mod(n_segments, 1) ~= 0 % Not an integer
    error('Something''s wrong, noninteger # of segments')
end

flip_flop = 1;
data1 = [];
data2 = [];
for i_segment = 1:n_segments
    curr = (i_segment - 1) * segment_length;
    % Prepare data, this is finding the correct index to segment
    inds = repmat({1}, 1, ndims(data));
    for i_dim = 1:ndims(data)
        inds{i_dim} = 1:size(data, i_dim);
    end
    
    working_dim = find(ismember(size(data), recording_length));
    inds{working_dim} = curr + 1: curr + segment_length;
    curr_data = data(inds{:});
    switch flip_flop
        case 1
            data1 = cat(working_dim, data1, curr_data);
            flip_flop = 2;
        case 2
            data2 = cat(working_dim, data2, curr_data);
            flip_flop = 1;
    end
end
end




