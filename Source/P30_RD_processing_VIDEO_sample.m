%% The following script is a sample script for processing video data acquired with a video camera

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
                     'video_camera',...
                     'WindTurbine');

path2video_list = {
                    fullfile(path2proj,'video_recording.mp4'), datetime(2022,09,14,12,39,18);...
                   };

azimuth2img = -21.95; % Degree
azimuzhOfReference = 36.47; % Degree
azimuth2turbine = azimuzhOfReference + azimuth2img; % degree

proc_para_list = {
                  'H115m', 1, azimuth2turbine;...
                  };

%% Processing
N_files = size(path2video_list,1);
for f_i = 1:N_files
    path2video = path2video_list{f_i,1};
    time_T0 = path2video_list{f_i,2};

    name_proc = proc_para_list{f_i,1};
    isIndirect = proc_para_list{f_i,2};
    azimuth = proc_para_list{f_i,3};

    video_tracking_v2(path2video,...
                      name_proc,...
                      isIndirect,...
                      azimuth,...
                      time_T0);
end