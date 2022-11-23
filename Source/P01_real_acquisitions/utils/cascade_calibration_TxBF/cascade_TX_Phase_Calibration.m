%  Copyright (C) 2020 Texas Instruments Incorporated - http://www.ti.com/
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
%
%

% cascade_TX_Phase_Calibration.m
%  
% Top level main test chain to perform TX phase shifter calibration. 
% Output is saved to calibrateTXPhaseResults.mat
% 
% Usage: modify the dataFolder_calib_data_path 


%% cascade_TX_Phase_Calibration.m
%
% The following script/function was modified from the above mentioned 
% TI-script. It can be used to process raw data acquired a TIDEP-01012 
% device. Specificially, a calibration dataset acquried with the 
% configuration of Cascade_Phase_Shifter_Calibration_AWRx.lua
% The result of this function is a MAT-file with phase shifts.
%
% -------------------------------------------------------------------------
%
% Input:
%   
% - path2proj:      The path to the main directory of the project. The
%                   directory needs to have a folder "01_Raw_Radar_Data" 
%                   with all the data acquired using the 
%                   Cascade_Phase_Shifter_Calibration_AWRx.lua 
%                   configuration.
%
% - trg_rng:        A value indicating the approximate distance of the
%                   corner cube to the radar sensor. [meter]
%
% - trg_rng_search: A value indicating the search distance (deviation from 
%                   inital trg_rng value). [meter]
%
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang, ETH Zürich (04nd August 2022)


function [calibrateTXPhaseResultsFile] = cascade_TX_Phase_Calibration(path2proj,trg_rng,trg_rng_search)

path2param = fullfile(path2proj,'00_Calibration_Files');
if ~isfolder(path2param)
    mkdir(path2param);
end

path2raw = fullfile(path2proj,'01_Raw_Radar_Data');

dataPlatform = 'TDA2';

DEBUG_PLOTS = 1;                % optionally display debug plots while calibration data is being processed, slow if activated!

numRX = 16;                     % number of MMWCAS-RF-EVM RX channels being processed, set to 16, for the full MMWCAS-RF-EVM RX channels
numPhaseShifterOffsets = 64;    % number of phase-shifter offset increments being processed, set to 64 for full phase-shifter range (number of datasets)
numChirpsLoopsPerFrame = 12;    % number of chirp-loops per frame. Only a single TX active per chirp-loop. 
numChirpsPerLoop = 64;          % number of chirps per TX active phase. 
searchBinsSkip = 0;             % number of bins to skip when looking for peak from corner reflector 

% select reference TX/Device channel for computing offsets - determined by 
% TX antenna array geometry and phase-shifter offset utilization (all RX
% channels are processed)
refTX = 1;
refPhaseOffset = 1;

fig(1) = figure('Name','Phase Offset',...
                'Units','centimeters',...
                'Position',[0,0,45,15]);  % plot phase offset values
ax(1) = axes;

fig(2) = figure('Name','Phase Offset Error',...
                'Units','centimeters',...
                'Position',[0,0,45,15]);  % plot phase offset error values
ax(2) = axes;

if DEBUG_PLOTS

    N_frames = numChirpsLoopsPerFrame*numRX*numPhaseShifterOffsets;

    fig(3) = figure('Name','Debugging',...
                    'Units','centimeters',...
                    'Position',[0,0,45,15],...
                    'visible','off'); 
    ax(3) = subplot(2,6,1:3); % FFT magnitude (dB)
    ax(4) = subplot(2,6,4:6); % FFT phase (deg)
    ax(5) = subplot(2,6,[7,8]); % accumulated target bin values
    ax(6) = subplot(2,6,[9,10]); % accumulated target distance values
    ax(7) = subplot(2,6,[11,12]); % accumulated target phase angle values

    %% Initialize video
    path4ani = fullfile(path2param,...
                        'calibrateTXPhaseResults_03.avi');
    VidRec = VideoWriter(path4ani);
    VidRec.FrameRate = 2;
    open(VidRec)

    idx_f = 1; % Index for Animation

    frame(N_frames) = getframe(fig(3));
end


%parameter file name for the test
pathGenParaFile = fullfile(path2param,'generateCalibrationMatrix_param.m');

%important to clear the same.m file, since Matlab does not clear cache
%automatically
clear(pathGenParaFile);

% loop through each folder in the calibration directory create 
% matrix with dimensions [devices, number of TX, number of phase shift offsets] 
% with the phase values for the detected range-bin peak
dataFolder_calib_data_info = dir(path2raw);

% find the first folder of cal data
for folderIdx = 1:length(dataFolder_calib_data_info)
    if(dataFolder_calib_data_info(folderIdx).isdir && dataFolder_calib_data_info(folderIdx).name ~= "." && dataFolder_calib_data_info(folderIdx).name ~= "..") 
        folderIdxStart = folderIdx;
        break;
    end

end

folderIdx = folderIdxStart; % start at first index that is a folder

phaseValuesBin = zeros(numChirpsLoopsPerFrame,numRX,numPhaseShifterOffsets);
phaseValuesTargetDistance = zeros(numChirpsLoopsPerFrame,numRX,numPhaseShifterOffsets);
phaseValues = zeros(numChirpsLoopsPerFrame,numRX,numPhaseShifterOffsets);

for idxPS = 1:numPhaseShifterOffsets % loop through phase shifter offset/transmitter/device

    if(dataFolder_calib_data_info(folderIdx).isdir) 
        disp(['dataFolder_calib_data_info(folderIdx).name = ', dataFolder_calib_data_info(folderIdx).name]);

        % create folder structure for this iteration        
        fileNameCascade.dataFolderName = strcat([dataFolder_calib_data_info(folderIdx).folder, '\', dataFolder_calib_data_info(folderIdx).name, '\']);
        [fileIdx_unique] = getUniqueFileIdx(fileNameCascade.dataFolderName);
        [fileNameStruct] = getBinFileNames_withIdx(fileNameCascade.dataFolderName, fileIdx_unique{1});                  

        %generate parameter file for the dataset
        parameter_file_gen_antennaCalib_json(fileNameCascade.dataFolderName, pathGenParaFile, dataPlatform);

        %generate calibration matrix object for the dataset
        genCalibrationMatrixObj = genCalibrationMatrixCascade(...
                    'pfile', pathGenParaFile,...
                    'calibrateFileName',fileNameCascade.dataFolderName,...
                    'targetRange', trg_rng);
        genCalibrationMatrixObj.binDataFile = fileNameStruct;
    end
            
    % use second frame for calibration 
    genCalibrationMatrixObj.frameIdx = 2;             

    % read in .bin files
    calData = cascade_Read_TX_Cal_Data(genCalibrationMatrixObj);
        
    if ~exist('trg_rng_search','var')
        searchRangeCCBin = 3;
        trg_rng_search = searchRangeCCBin * genCalibrationMatrixObj.rangeResolution; % distance in [m] for searching the peak from the corner reflector
    end

    searchRangeCCBin = ceil(trg_rng_search/genCalibrationMatrixObj.rangeResolution);
    searchCCBin = round(trg_rng/genCalibrationMatrixObj.rangeResolution);

    for idxTX = 1:numChirpsLoopsPerFrame % loop through each TX phase   
        for idxRX = 1:numRX % loop through each RX
        
            calData_1DFFT(:, :, idxTX, idxRX) = fftshift(fft(calData(:, :, idxRX, idxTX), genCalibrationMatrixObj.numSamplePerChirp, 1));
            calData_2DFFT(:, :, idxTX, idxRX) = 1/(numChirpsPerLoop) * fftshift(fft(calData_1DFFT(:, :, idxTX, idxRX), numChirpsPerLoop, 2));
            
          
            % find target peak bin in 2D-FFT 0-velocity bin (skip close
            % bins to avoid DC leakage or bumper reflections

            calData_1DFFT_tmp = abs(squeeze(calData_2DFFT(:, numChirpsPerLoop/2 + 1, idxTX, idxRX)));
            TargetBinValue = max(calData_1DFFT_tmp(searchCCBin-searchRangeCCBin:searchCCBin+searchRangeCCBin));
            TargetBinIdx = find(calData_1DFFT_tmp==TargetBinValue);
            phaseValuesBin(idxTX, idxRX, idxPS) = TargetBinIdx + searchBinsSkip - 1;
            phaseValuesTargetDistance(idxTX, idxRX, idxPS) = TargetBinIdx * genCalibrationMatrixObj.rangeResolution;

            % record phase at target peak bin
            phaseValues(idxTX, idxRX, idxPS) = angle(squeeze(calData_2DFFT(phaseValuesBin(idxTX, idxRX, idxPS), numChirpsPerLoop/2 + 1, idxTX, idxRX))) * 180 / pi; 

            % debug plots 
            if DEBUG_PLOTS
                numBins = size(calData_2DFFT,1);

                plot(ax(3), 10*log(abs(squeeze(calData_2DFFT(:, numChirpsPerLoop/2 + 1, idxTX, idxRX)))));
                hold(ax(3), 'on');
                plot(ax(3), phaseValuesBin(idxTX, idxRX, idxPS), 10*log(abs(TargetBinValue)), '-o', 'color', 'red');
                hold(ax(3), 'off');
                title(ax(3), 'Calibration Target Power vs. IF bins');
                xlabel(ax(3), '1D-FFT Spectrum (bins)'); 
                ylabel(ax(3), '1D-FFT Magnitude (dB)');
                ylim(ax(3),[0,150]);
                xlim(ax(3),[0,numBins]);
                hold on
                xline(ax(3),trg_rng/genCalibrationMatrixObj.rangeResolution);
                hold off

                plot(ax(4), angle(squeeze(calData_2DFFT(:, numChirpsPerLoop/2 + 1, idxTX, idxRX))) * 180 / pi );
                hold(ax(4), 'on');
                plot(ax(4), phaseValuesBin(idxTX, idxRX, idxPS), phaseValues(idxTX, idxRX, idxPS) , '-o', 'color', 'red');
                hold(ax(4), 'off');                
                title(ax(4), 'Calibration Target Phase vs. IF bins');
                xlabel(ax(4), '1D-FFT Spectrum (bins)');
                ylabel(ax(4), '1D-FFT Phase (degrees)');
                ylim(ax(4),[-180,180]);
                xlim(ax(4),[0,numBins]);

                scatter(ax(5), ...
                        1:numPhaseShifterOffsets,...
                        squeeze(phaseValuesBin(idxTX, idxRX, :)),...
                        'filled');                
                title(ax(5), 'Calibration Target Detected Index');
                xlabel(ax(5), 'Phase Shifter Offset (5.625 degrees/LSB)');
                ylabel(ax(5), 'Calibration Target Sampled IF Index');
                hold on
                yline(ax(5),trg_rng/genCalibrationMatrixObj.rangeResolution);
                hold off
                ylim(ax(5),[0,numBins]);
                xlim(ax(5),[0,numPhaseShifterOffsets]);

                scatter(ax(6), ...
                        1:numPhaseShifterOffsets,...
                        squeeze(phaseValuesTargetDistance(idxTX, idxRX, :)),...
                        'filled');
                title(ax(6), 'Calibration Target Distance');
                xlabel(ax(6), 'Phase Shifter Offset (5.625 degrees/LSB)');
                ylabel(ax(6), 'Target Distance (meters)');
                ylim(ax(6),[0,genCalibrationMatrixObj.rangeResolution*numBins]);
                xlim(ax(6),[0,numPhaseShifterOffsets]);
                hold on
                yline(ax(6),trg_rng);
                hold off

                plot(ax(7), squeeze(phaseValues(idxTX, idxRX, :)));
                title(ax(7), 'Calibration Target Phase vs. IF bins');
                xlabel(ax(7), 'Phase Shifter Offset (5.625 degrees/LSB)');
                ylabel(ax(7), 'Target 1D-FFT Phase (degrees)');            
                ylim(ax(7),[-180,180]);


                sgtitle(fig(3),sprintf('Phase Shifter % 2d of TX%02d/RX%02d',...
                                idxPS,...
                                idxTX,...
                                idxRX));

                frame(idx_f) = getframe(fig(3));
                
                if idx_f == 1
                    wbar = waitbar(idx_f/N_frames,'Processing Visualisation...');
                else
                    waitbar(idx_f/N_frames,wbar,sprintf('Processing Visualisation... (Frame: % 5d / %d)',idx_f,N_frames));
                end

                idx_f = idx_f + 1;

            end
            
        end   

    end

    folderIdx = folderIdx + 1;

end


if DEBUG_PLOTS
    writeVideo(VidRec, frame);
    close(VidRec);
    close(wbar);
end

% compute phase offsets and phase errors
for idxTX = 1:numChirpsLoopsPerFrame % loop through transmitter/device
    for idxRX = 1:numRX % loop through transmitter/device
        for idxPS = 1:numPhaseShifterOffsets % loop through phase shifter offset/transmitter/device

            % reference the phase-shifter 0 setting, from channel refDevice, refTX 
            % as reference and compute offset phase values for the other settings
            phaseOffsetValues(idxTX, idxRX, idxPS) = phaseValues(idxTX, idxRX, idxPS) - phaseValues(refTX, idxRX, refPhaseOffset); 
            % find error between actual phase offset and expected (ideal) phase offset
            % this error will be a combination of static phase error and phase-shifter induced error
            phaseOffsetError(idxTX, idxRX, idxPS) = phaseValues(idxTX, idxRX, idxPS) - idxPS * 5.625;

        end
    end
end




 
% phase_shift_values = 0:5.625:5.625*64-1;

for TXIdx = 1:numChirpsLoopsPerFrame % loop through transmitter/device
    for RXIdx = 2:2 % loop through receivers   

        plot(ax(1), squeeze(phaseValues(TXIdx, RXIdx, :)), 'color',[TXIdx/numChirpsLoopsPerFrame, TXIdx/numChirpsLoopsPerFrame,  RXIdx * 0.0625]);
        hold(ax(1), 'on');

    end
end
       
%% plot results
title(ax(1), 'TX Channel Target Phase Value (deg) vs. TX Phase-Shifter Offset Value (5.625 deg increment)');
xlabel(ax(1), 'TX Phase-Shifter Offset Value (5.625 deg increment)'); 
ylabel(ax(1), 'TX Channel Absolute Phase (deg)'); 
legend(ax(1), "TX1", "TX2", "TX3", "TX4", "TX5", "TX6", "TX7", "TX8", "TX9", "TX10", "TX11", "TX12");
grid(ax(1), 'on');

for TXIdx = 1:numChirpsLoopsPerFrame % loop through transmitter/device
    for RXIdx = 1:1 % loop through receivers    

        plot(ax(2), squeeze(phaseOffsetValues(TXIdx, RXIdx, :)),...
                'Color',[TXIdx/numChirpsLoopsPerFrame, TXIdx/numChirpsLoopsPerFrame,  RXIdx * 0.0625]);
        hold(ax(2), 'on');

    end
end

title(ax(2), 'TX Channel Phase Offset (deg) vs. TX Phase-Shifter Offset Value (5.625 deg increment)');
xlabel(ax(2), 'TX Phase-Shifter Offset Value (5.625 deg increment)'); 
ylabel(ax(2), 'TX Channel Phase Offset (deg)'); 
legend(ax(2), "TX1", "TX2", "TX3", "TX4", "TX5", "TX6", "TX7", "TX8", "TX9", "TX10", "TX11", "TX12");
grid(ax(2), 'on');
    

%% save results
% save calibration data to .mat file
calibrateTXPhaseResultsFile = fullfile(path2param,...
                                       'calibrateTXPhaseResults.mat');

save(calibrateTXPhaseResultsFile, ...
    'phaseOffsetValues',...
    'phaseValues',...
    'phaseOffsetError');

for fi = 1:length(fig)
    path4fig = fullfile(path2param,...
                        sprintf('calibrateTXPhaseResults_%02d.fig',fi));
    
    savefig(fig(fi),path4fig);
    exportgraphics(fig(fi),...
                    replace(path4fig,'.fig','.png'),...
                    'Resolution',600);
end


end


