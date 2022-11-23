%% The following script is a sample script for processing radar data acquired with a TIDEP-01012 system

%% Clear everything
close all
clearvars
clc

%% Initalisation
[file_dir,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(fullfile(file_dir)));

% Define Project
proj_list = {'MIMO_A77_Calandawind_20220914_124115_Visualisation';}; % Project name consistent with <project_name>

filt_by_dist_list = {{1,2,200};}; % Distance Limitation {True/False, D_Min, D_Max}

filt_by_azi_list = {{1,-65,65};}; % Azimuth Limitation {True/False, Az_Min, Az_Max}

filt_by_asi_list = {{1,-inf};}; % Filter by Amplitude Stability Index {True/False, ASI_Min}

filt_by_lr_list = {{1,-500,500};}; % Cross Range Limitation {True/False, CR_Min, CR_Max}

filt_by_coh_list = {{1,-inf};}; % Filter by Coherence {True/False, COH_Min}

filt_by_maxDisp_list = {{1,-inf,inf};}; % Filter by Maximum Displacement {True/False, dDNeg_max, dDPos_max}

filt_by_aoi_list = {{1,1,[0]};}; % Filter by Area of Interest {True/False, Number of AoI, isCircle True/False}

time_select_list = {{1, datetime(2022,09,14,11,12,45),datetime(2022,09,14,11,13,15)};}; % Filter by Time {True/False, T_Start, T_End}

time_select_zoom_list = {{1, datetime(2022,09,14,11,12,45),datetime(2022,09,14,11,13,15)};}; % Filter by Time (Zoom for psi2timeseries function)

filt_by_aoi_zoom_list = {{1,1,[0]};}; % Filter by Area of Interest (Zoom for psi2timeseries function)

for p_i = 1:length(proj_list)

    name2proj = proj_list{p_i};
    path2proj = fullfile(file_dir,...
                         'D00_sample_data',...
                         'real',...
                         name2proj);
    
    %% Processing: From Raw Data to SLC
    filt_by_rng  = filt_by_dist_list{p_i};
    filt_by_azi  = filt_by_azi_list{p_i};
    filt_by_asi  = filt_by_asi_list{p_i};
    
    cascade_MIMO_01_raw2slc(path2proj,...
                            filt_by_rng,...
                            filt_by_azi,...
                            filt_by_asi);

    %% Processing: From SLC Data to PSI
    % Settings for Geometrical Filtering
    filt_by_aoi = filt_by_aoi_list{p_i};
    filt_by_rng = filt_by_dist_list{p_i};
    filt_by_lr  = filt_by_lr_list{p_i};
    filt_by_azi = filt_by_azi_list{p_i};
    
    % Settings for Statistical Filtering
    filt_by_asi = filt_by_asi_list{p_i};
    filt_by_coh = filt_by_coh_list{p_i};
    filt_by_maxDisp = filt_by_maxDisp_list{p_i};
    
    % Settings for Temporal Filtering
    filt_by_time = time_select_list{p_i};
    
    cascade_MIMO_02_slc2psi(path2proj,...
                           filt_by_rng,...
                           filt_by_lr,...
                           filt_by_azi,...
                           filt_by_asi,...
                           filt_by_coh,...
                           filt_by_maxDisp,...
                           filt_by_time,...
                           filt_by_aoi);


    %% Processing: From PSI to Coordinate Components
    % Settings for Temporal Filtering
    filt_by_asi = filt_by_asi_list{p_i};
    filt_by_time = time_select_zoom_list{p_i};
    create_aoi = filt_by_aoi_zoom_list{p_i}{2};

    cascade_MIMO_03_psi2timeseries(path2proj, ...
                                    filt_by_time, ...
                                    filt_by_asi, ...
                                    create_aoi);

end