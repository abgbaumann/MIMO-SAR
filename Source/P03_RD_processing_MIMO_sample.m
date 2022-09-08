%% The following script is a sample script for processing radar data 
%  acquired with a TIDEP-01012 system

%% Clear everything
close all
clearvars
clc

%% Initalisation
[file_dir,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(fullfile(file_dir)));

% Define Project
proj_list = {'MIMO_A77_Calandawind_20220817_125659_00120000ms';...
             'MIMO_A77_Calandawind_20220817_131729_00120000ms';...
             'MIMO_A77_Calandawind_20220817_134511_00120000ms';...
             'MIMO_A77_Calandawind_20220817_145207_00120000ms';...
             'MIMO_B77_Calandawind_20220817_125735_00120000ms';...
             'MIMO_B77_Calandawind_20220817_131752_00120000ms';...
             'MIMO_B77_Calandawind_20220817_134555_00120000ms';...
             'MIMO_C77_Calandawind_20220817_130109_00120000ms';...
             'MIMO_C77_Calandawind_20220817_131859_00120000ms';...
             'MIMO_C77_Calandawind_20220817_134713_00120000ms';...
             };


filt_by_dist_list = {{1,130,220};...
                     {1,130,220};...
                     {1,100,175};...
                     {1,60,160};...
                     {1,108,225};...
                     {1,108,225};...
                     {1,108,225};...
                     {1,50,190};...
                     {1,50,190};...
                     {1,50,190};...
                     }; % Calandawind A

filt_by_azi_list = {{1,-30,30};...
                    {1,-30,30};...
                    {1,-60,60};...
                    {1,-30,30};...
                    {1,-26,34};...
                    {1,-26,34};...
                    {1,-26,34};...
                    {1,-28,28};...
                    {1,-28,28};...
                    {1,-28,28};...
                    }; % Calandawind A

filt_by_asi_list = {{0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    }; % Calandawind A

filt_by_lr_list = {{0,-100,100};...
                   {0,-100,100};...
                   {0,-100,100};...
                   {0,-80,80};...
                   {0,-100,100};...
                   {0,-100,100};...
                   {0,-100,100};...
                   {0,-80,80};...
                   {0,-80,80};...
                   {0,-80,80};...
                   }; % Calandawind A

filt_by_coh_list = {{0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    {0,-inf};...
                    };

filt_by_maxDisp_list = {{0,-inf,inf};...
                        {0,-inf,inf};...
                        {0,-inf,inf};...
                        {0,-inf,inf};...
                        {0,-inf,inf};...
                        {0,-inf,inf};...
                        {0,-inf,inf};...
                        {0,-inf,inf};...
                        {0,-inf,inf};...
                        {0,-inf,inf};...
                        };

filt_by_aoi_list = {{1,1,[0]};...
                    {1,1,[0]};...
                    {1,1,[0]};...
                    {1,1,[1]};...
                    {1,1,[0]};...
                    {1,1,[0]};...
                    {1,1,[0]};...
                    {1,1,[1]};...
                    {1,1,[1]};...
                    {1,1,[1]};...
                    };

time_select_list = {{1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                    {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                    {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                    {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                    {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                    {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                    {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                    {1, datetime(2022,08,17,13,47,00,000),datetime(2022,08,17,13,50,00,000)};...
                    {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                    {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                    };

time_select_zoom_list = {{1, datetime(2022,08,17,13,01,00,000),datetime(2022,08,17,13,03,00,000)};...
                         {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                         {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                         {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                         {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                         {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                         {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                         {1, datetime(2022,08,17,13,47,00,000),datetime(2022,08,17,13,50,00,000)};...
                         {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                         {1, datetime(2022,08,17,00,00,00,000),datetime(2022,08,17,23,59,59,000)};...
                        };

for p_i = 8%:1:length(proj_list)%[1]%[1,5,8,2,6,9,3,7,10,4]%1:1:length(proj_list)

    name2proj = proj_list{p_i};
    path2proj = fullfile(file_dir,'D00_sample_data','real',name2proj);
    
    %% Processing: From Raw Data to SLC
%     filt_by_rng  = filt_by_dist_list{p_i};  % Filter by Distance: {true/false, minD, maxD}
%     filt_by_azi  = filt_by_azi_list{p_i};   % Filter by Azimuth: {true/false, minAz, maxAz}
%     filt_by_asi  = filt_by_asi_list{p_i};   % Filter by Amplitude Stability Index: [0,1];
%     
%     cascade_MIMO_01_raw2slc(path2proj,filt_by_rng,filt_by_azi,filt_by_asi);
%     
%     close all

    %% Processing: From SLC Data to PSI
    % Settings for Geometrical Filtering (set first value of each cell to 0 to
    % deactivate the filtering)
    filt_by_aoi = filt_by_aoi_list{p_i};    % Keep everything within the Area of Interest [-]
    filt_by_rng = filt_by_dist_list{p_i};   % Keep everything within the Distance ranges [m]
    filt_by_lr  = filt_by_lr_list{p_i};     % Keep everything within the Cross Range ranges [m]
    filt_by_azi = filt_by_azi_list{p_i};    % Keep everything within the Azimuth ranges [deg]
    
    % Settings for Statistical Filtering (set first value of each cell to 0 to
    % deactivate the filtering)
    filt_by_asi = filt_by_asi_list{p_i};   % Keep everything with a higher Amplitude Stability Index [-]
    filt_by_coh = filt_by_coh_list{p_i};    % Keep everything with a higher Coherence [-] as derived from the 3x3 Neighbourhood
    filt_by_maxDisp = filt_by_maxDisp_list{p_i}; % Keep everything with a smaller displacement
    
    % Settings for Temporal Filtering (set first value of each cell to 0 to
    % deactivate the filtering)
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
    % Settings for Temporal Filtering (set first value of each cell to 0 to
    % deactivate the filtering)
%     filt_by_asi = filt_by_asi_list{p_i};
%     filt_by_time = time_select_zoom_list{p_i};
%     create_aoi = 4;
% 
%     cascade_MIMO_03_psi2timeseries(path2proj,filt_by_time,filt_by_asi,create_aoi);

end