%% The following script is a sample script for processing radar data 
%  acquired with a TIDEP-01012 system to calibrate the sensor for
%  using in the TX beamforming mode (TXBF).

%% Clear everything
close all
clearvars
clc

%% Initalisation
[file_dir,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(fullfile(file_dir)));

proj_list = {'TxBF_A77_phaseShiftCalibration_20220811';...
             'TxBF_B77_phaseShiftCalibration_20220811'};

% trg_rng_list = [113.5,...
%                 113.5]; % [m]
trg_rng_list = [44.5,...
                44.5]; % [m]

%% Processing Beamforming Calibration
for p_i = 1:length(proj_list)

    proj_name = proj_list{p_i};
    path2proj = fullfile(file_dir,'D00_sample_data','real',proj_name);

    trg_rng = trg_rng_list(p_i); % estimated corner-reflector target range (in meters) 
    
    path2cal = cascade_TX_Phase_Calibration(path2proj, trg_rng);
    
    [path2cal_1,path2cal_2] = TXBF_PS_LUT_Generate(path2cal);
    
    path2cal_3 = TXBF_Create_PSCal_Advanced_Frame_Config(path2cal); % Matrix for TxBF - Acquisition
    
    close all

end

