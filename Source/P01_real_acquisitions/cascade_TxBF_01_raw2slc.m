%  Copyright (C) 2018 Texas Instruments Incorporated - http://www.ti.com/
%
%
%   Redistribution and use in source and binary forms, with or without
%   modification, are permitted provided that the following conditions
%   are met:
%
%     Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer.
%
%     Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the
%     distribution.
%
%     Neither the name of Texas Instruments Incorporated nor the names of
%     its contributors may be used to endorse or promote products derived
%     from this software without specific prior written permission.
%
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
%   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
%   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
%   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
%   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
%   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% read the adc data with advanced frame configuration with TX beamforming,
% stich the multiple scans together with RX beamforming

%read the adc data and postprocess it
%path_in = 'Z:\Research\04_InternalPhD\Baumann_Andreas\50_Processing_Final\D00_sample_data\real\Cascade_TxBF_A77_ETH_Buildings\01_Raw_Radar_Data\';
%path2calib = 'Z:\Research\04_InternalPhD\Baumann_Andreas\50_Processing_Final\D00_sample_data\real\MIMO_A77_phaseMismatchCalibration_LR230\00_Calibration_Files\phaseMismatchCalibration.mat';
% angle2steer = [-45:2:45]


function [] = cascade_TxBF_01_raw2slc(path2proj, path2calib, angle2steer, filt_by_dist, filt_by_azi, filt_by_asi)

[~,proj_name]= fileparts(fileparts(path2proj));
path_in = fullfile(path2proj,'01_Raw_Radar_Data'); % Input Path to Raw Data
path_out = fullfile(path2proj,'02_SLC_Radar_Data'); % Output Path to SLC Data

if ~isfolder(path_out)
    mkdir(path_out);
end

 %% Start Processing
step_i = 1; % Numbering for command window message
fprintf('Step %02d: Process Radar Data (%s)\n', step_i, proj_name);

%% define what angles to steer in TX beamforming mode
paramsConfig.anglesToSteer = angle2steer; 

% Only Frame based supported in the Post-processing script
paramsConfig.Chirp_Frame_BF = 0; % 1 - chirp based beam steering, 0 - frame based beam steering

paramsConfig.NumAnglesToSweep = length(paramsConfig.anglesToSteer);
numAngles = paramsConfig.NumAnglesToSweep;

%parse Json file
paramsConfig = parameter_gen_from_Jason(path_in, paramsConfig);

%load MIMO calbibration params
load(path2calib);
paramsConfig.Slope_MHzperus_calib = params.Slope_MHzperus;
BF_MIMO_ref = calibResult.RxMismatch;

%re-arrange calibration data
RangeMat = calibResult.RangeMat;
RangeMat = RangeMat(:, paramsConfig.TI_Cascade_RX_ID);
PeakValMat = calibResult.PeakValMat;
PeakValMat = PeakValMat(:, paramsConfig.TI_Cascade_RX_ID);

%pass parameters from paramsConfig
numFrames = paramsConfig.Num_Frames;
Rx_Ant_Arr = paramsConfig.TI_Cascade_RX_ID;
anglesToSteer = paramsConfig.anglesToSteer;

currentFolder = pwd;
cd(path_in);
listing = dir('*_data.bin');
cd(currentFolder);

% Get Unique File Idxs in the "dataFolder_test"   
[fileIdx_unique] = getUniqueFileIdx(path_in);
N_file_sets = length(fileIdx_unique);

frameCountGlobal    = 0;
t0 = tic;
totNumFrames = 0;
for i_file = 1:N_file_sets
    [fileNameStruct]= getBinFileNames_withIdx(path_in, fileIdx_unique{i_file});
    [totNumFramesTmp, ~] = getValidNumFrames(fullfile(path_in, fileNameStruct.masterIdxFile));
    totNumFrames = totNumFrames + totNumFramesTmp;
    if i_file == 1
        totNumFrames = totNumFrames - 1;
    end
end
for i_file = 1:N_file_sets
    
    % Get File Names for the Master, Slave1, Slave2, Slave3   
    [fileNameStruct]= getBinFileNames_withIdx(path_in, fileIdx_unique{i_file});
    
    RX_TXBF_calVec = ones(1,length(paramsConfig.Rx_Elements_To_Capture));
    
    % Get Valid Number of Frames 
    [numValidFrames dataFileSize] = getValidNumFrames(fullfile(path_in, fileNameStruct.masterIdxFile));
    
    % intentionally skip the first frame due to TDA2 
    if i_file == 1 % Skip the first frame of the dataset
        startValidFrame = 2;
    else
        startValidFrame = 1;
    end

    no_of_local_acquisitions = numValidFrames-(startValidFrame-1);

    step_i = step_i + 1;

    for frameIdx = startValidFrame:numValidFrames

        frameCountGlobal = frameCountGlobal+1;
        frameCountLocal = frameIdx-(startValidFrame-1);

        paramsConfig.frameId = frameIdx;
        [radar_data_TXBF]= fcn_read_AdvFrmConfig_BF_Json(fileNameStruct,paramsConfig);
        
        % TI_Cascade_RX_ID: RX channel ID order to be read in sequence
        radar_data_TXBF = radar_data_TXBF(:,:,paramsConfig.TI_Cascade_RX_ID,:);
        
        %% ------------------process TXBF data
        startRangeInd = 5;
        lastRangeIndToThrow = 20;
        
        [ complex_data_static_tmp, y_axis, x_axis] = Plot_advFraConfig_TXBF_rangeAzimuth_stich(radar_data_TXBF,...
            Rx_Ant_Arr,paramsConfig,anglesToSteer,BF_MIMO_ref(Rx_Ant_Arr), startRangeInd, lastRangeIndToThrow);

        if frameCountLocal == 1
            complex_data_static = zeros([size(complex_data_static_tmp),no_of_local_acquisitions]);
        end
        
        complex_data_static(:,:,frameCountLocal) = complex_data_static_tmp;

        time_i = toc(t0)/frameCountGlobal;
        time_rem = time_i*(totNumFrames-frameCountGlobal);
        time_rem_min = floor(time_rem/60);
        time_rem_sec = floor(time_rem - (time_rem_min*60));
        time_end = datetime() + seconds(time_rem);
    
        if exist('lineLength','var')
            fprintf(repmat('\b',1,lineLength));
        end
    
        print_str = sprintf("Step %02d: Frame % 5d processed (SLC).",...
                            step_i,...
                            frameCountGlobal);
        if frameCountGlobal/totNumFrames<0.1
            lineLength = fprintf("%s\n",print_str);
        else
            lineLength = fprintf("%s (Ending in %d:%02dmin at %s)\n",...
                    print_str,...
                    time_rem_min,...
                    time_rem_sec,...
                    datestr(time_end,"dd.mm.yyyy HH:MM:SS"));
        end
    end    

    clearvars lineLength
    
    ind = strfind(path_in, filesep);
    testName = path_in(ind(end-1)+1:(ind(end)-1));
    path_mat_03 = fullfile(path_out,sprintf("%s_%05d.mat",testName,i_file)); 
    
    %% Calculate Coherence of first to second, to N/2, and to N
    step_i = step_i + 1;
    fprintf('Step %02d: Calculate Coherences.\n',step_i);
        
    step_i = step_i + 1;
    coherence = zeros(size(complex_data_static));
    N_slc = size(complex_data_static,3);
    for slc_i = 1:N_slc
        coherence(:,:,slc_i) = get_coherence(complex_data_static(:,:,1),...
                                        complex_data_static(:,:,slc_i),...
                                        3);
    
        if exist('lineLength','var')
            fprintf(repmat('\b',1,lineLength));
        end
    
        print_str = sprintf("Step %02d: Frame % 5d / %d processed (Coherence).",...
                            step_i,...
                            slc_i,...
                            N_slc);
        lineLength = fprintf("%s\n",print_str);
    end
    
    coh_mean = mean(coherence,3);
    coh_std = std(coherence,[],3);
    
    coh_class_tresh = 0.9;
    coh_mat_e_filt = ones(size(coh_mean));
    coh_mat_e_filt(squeeze(coherence(:,:,end))<=coh_class_tresh) = 0;
    coh_mat_s_filt = 2 * ones(size(coh_mean));
    coh_mat_s_filt(squeeze(coherence(:,:,2))<=coh_class_tresh) = 0;
    coh_class = coh_mat_e_filt-coh_mat_s_filt;
    
    %% Visualise Data Amplitudes and Coherences
    step_i = step_i + 1;
    fprintf('Step %02d: Visualise Amplitudes and Coherences\n',step_i);
    close all
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    amp = abs(complex_data_static);
    
    % Amplitude Image
    ax(1)=subplot(2,2,1);
    plot_polar_range_azimuth_2D_AB_preAX(y_axis,...
                                            x_axis,...
                                            log10(mean(amp,3)));
    hcb = colorbar;
    hcb.Label.String = "Amplitude [log10(A)]";
    
    ax(2)=subplot(2,2,2);
    dim=ndims(complex_data_static);
    n_max = 4000;
    if size(complex_data_static,dim)>n_max
        n_skip = ceil(size(complex_data_static,dim)/n_max);
    else
        n_skip = 1;
    end
    if dim==3
        amptmp = squeeze(amp(:,:,1:n_skip:end));
    else
        error('Data format for Amplitude with given Dimension is unknown!\n');
    end
    amp_sigma = std(amptmp,[],dim);
    amp_mu = mean(amptmp,dim);
    plt_asi = 1 - (amp_sigma ./ amp_mu);
    plt_asi(plt_asi<=0) = 0;
    plot_polar_range_azimuth_2D_AB_preAX(y_axis,...
                                          x_axis,...
                                          plt_asi);
    hcb = colorbar;
    hcb.Label.String = "Amplitude Stability Index [-]";
    climits = [prctile(plt_asi(:),50),prctile(plt_asi(:),95)];
    clim(climits);
    caxis(climits);
    
    ax(3)=subplot(2,2,3);
    plot_polar_range_azimuth_2D_AB_preAX(y_axis,...
                                          x_axis,...
                                          coh_mean);
    hcb = colorbar;
    hcb.Label.String = "Mean Coherence [0...1]";
    climits = [prctile(coh_mean(:),50),prctile(coh_mean(:),90)];
    clim(climits);
    caxis(climits);
    
    ax(4)=subplot(2,2,4);
    plot_polar_range_azimuth_2D_AB_preAX(y_axis,...
                                          x_axis,...
                                          coh_std);
    colormap(ax(4), flipud(colormap))
    hcb = colorbar;
    hcb.Label.String = "Standard Deviation of Coherence";
    climits = [prctile(coh_std(:),5),prctile(coh_std(:),50)];
    clim(climits);
    caxis(climits);
    
    Link = linkprop(ax,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim'});
    setappdata(gcf, 'StoreTheLink', Link);
    
    sgtitle(sprintf('%s - Amplitude and Coherence (%d)',replace(proj_name,'_',' '),i_file));
    
    exportgraphics(fig,replace(path_mat_03,'.mat','_1.png'),'Resolution',600);
    savefig(fig,replace(path_mat_03,'.mat','_1.fig'));

    %% Saving Data to MAT-File.
    step_i = step_i + 1;
    fprintf('Step %02d: Compress Data by Filtering based:\n',step_i);
        
    if filt_by_dist{1}
        fprintf('- Range Limits (%.2f - %.2f m)\n',filt_by_dist{2},filt_by_dist{3});
        rng = sqrt(x_axis(1,:).^2+y_axis(1,:).^2);
        idx_rmv = rng<filt_by_dist{2} | rng>filt_by_dist{3};
        x_axis(:,idx_rmv) = []; 
        y_axis(:,idx_rmv) = []; 
        complex_data_static(:,idx_rmv,:) = [];
        coherence(:,idx_rmv,:) = [];
        coh_mean(:,idx_rmv) = [];
        coh_std(:,idx_rmv) = [];
        coh_class(:,idx_rmv) = [];
    end
        
    if filt_by_azi{1}
        fprintf('- Azimuth Limits (%.1f - %.1f deg)\n',filt_by_azi{2},filt_by_azi{3});
        azi = atan2(y_axis(:,1),x_axis(:,1)) * 180 / pi;
        idx_rmv = azi<filt_by_azi{2} | azi>filt_by_azi{3};
        x_axis(idx_rmv,:) = []; 
        y_axis(idx_rmv,:) = []; 
        complex_data_static(idx_rmv,:,:) = [];
        coherence(idx_rmv,:,:) = [];
        coh_mean(idx_rmv,:) = [];
        coh_std(idx_rmv,:) = [];
        coh_class(idx_rmv,:) = [];
    end
        
        x_axis_coh = x_axis;
        y_axis_coh = y_axis;
        
    if filt_by_asi{1} == 1
        
        if length(filt_by_asi{2}) == 1
            fprintf('- Amplitude Stability Index (>=%.2f)\n',filt_by_asi{2});
        end
        
         [complex_data_static,x_axis,y_axis,filt_by_asi{2},idx_a,idx_b] = ...
            filter_by_asi(complex_data_static,...
                          x_axis,...
                          y_axis,...
                          i_file,...
                          filt_by_asi{2});
        
        % Update Coherences
        coherence_tmp = zeros(size(complex_data_static));
        coh_mean_tmp = zeros(length(idx_a),1);
        coh_std_tmp = coh_mean_tmp;
        coh_class_tmp = coh_mean_tmp;
        for idx_i = 1:length(idx_a)
            coherence_tmp(idx_i,:) = coherence(idx_a(idx_i),idx_b(idx_i),:);
            coh_mean_tmp(idx_i,:) = coh_mean(idx_a(idx_i),idx_b(idx_i));
            coh_std_tmp(idx_i,:) = coh_std(idx_a(idx_i),idx_b(idx_i));
            coh_class_tmp(idx_i,:) = coh_class(idx_a(idx_i),idx_b(idx_i));
        end
        
        coherence = coherence_tmp;
        coh_mean = coh_mean_tmp;
        coh_std = coh_std_tmp;
        coh_class = coh_class_tmp;
        x_axis_coh = x_axis;
        y_axis_coh = y_axis;
    end
    
    step_i = step_i + 1;
    fprintf('Step %02d: Visualise Filtered Image: Amplitude\n',step_i);
    
    fig2 = figure('units','normalized','outerposition',[0 0 1 1]);
    amp = abs(complex_data_static);
    r=sqrt(y_axis.^2+x_axis.^2)*1.4*pi/180*100;
    scatter(y_axis,x_axis,r,log10(mean(amp,2)),'filled');
    hcb = colorbar;
    hcb.Label.String = "Amplitude [log10(A)]";
    axis equal
    box on
    xlabel('Cross Range [m]');
    ylabel('Along Range [m]');
    sgtitle(sprintf('%s - Filtered Amplitude Image (%d)',replace(proj_name,'_',' '),i_file));
    
    exportgraphics(fig2,replace(path_mat_03,'.mat','_2.png'),'Resolution',600);
    savefig(fig2,replace(path_mat_03,'.mat','_2.fig'));
    
    coh.class = coh_class;
    coh.class_id = [-2,-1,0,1];
    coh.class_descr = {'Higher Coherence at Start',...
                             'High Coherence',...
                             'Low Coherence',...
                             'Higher Coherence at End'};
    
    coh.coh_1N = coherence;
    coh.mean   = coh_mean;
    coh.std    = coh_std;
    coh.x_axis = x_axis_coh;
    coh.y_axis = y_axis_coh;
    
    step_i = step_i + 1;
    fprintf('Step %02d: Saving Data to MAT-File\n',step_i);
    save(path_mat_03,...
        'complex_data_static',...
        'coh',...
        'x_axis',...
        'y_axis',...
        'numAngles',...
        '-v7.3');

    clearvars lineLength

end

%%
end