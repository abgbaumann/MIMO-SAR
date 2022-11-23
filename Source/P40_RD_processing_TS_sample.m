%% The following script is a sample script for processing data acquired by a total station

%% Clear everything
close all
clearvars
clc

%% Initalisation
[file_dir,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(fullfile(file_dir)));

%% Processing Parameters
path2proj = fullfile(file_dir,...
                     'D00_sample_data',...
                     'total_station',...
                     'WindTurbine');
input_ts_list = {
                fullfile(path2proj,'total_station_raw.txt'), 2022, -duration(0,21,38,936), 172.3747;...
               }; % Input Information: path to raw data, year of acquisition, time shift in case of time error, azimuth of acquistion


%% Processing
N_files = size(input_ts_list,1);
for f_i = 1:N_files
    path2ts = input_ts_list{f_i,1};
    year = input_ts_list{f_i,2};
    time_shift =  input_ts_list{f_i,3};
    azimuth = input_ts_list{f_i,4};

    data_ts = total_station_tracking(path2ts, year, time_shift, azimuth);    
end