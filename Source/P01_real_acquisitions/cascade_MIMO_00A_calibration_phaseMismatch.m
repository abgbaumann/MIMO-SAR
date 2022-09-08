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
%
%

% cascade_MIMO_antennaCalib.m
%  
% Top level main test chain to perform antenna calibration. 
% genCalibrationMatrixCascade module is first initialized before
% actually used in the chain. The output is saved to a .mat file to be used
% later in data processing

%% cascade_MIMO_00A_calbration_antenna.m
%
% The following scrip/function was modified from the above mentioned 
% TI-scirpt. It can be used to process raw data acquired with a TIDEP-01012 
% device. The result of this function is a dataset consisting of SLCs.
%
% -------------------------------------------------------------------------
%
% Input:
%   
% - tbd:            tbd.
%                   The directory should have two folders.
%                   I) In "00_Parameter_Files":
%                    a) A sensor calibration file (Sensor_Calibration_at_77GHz.mat)
%                    b) A parameter file (module_param_1Chirp.m)
%                   II) In "01_Raw_Radar_Data":
%                    a) All acquisition files exported from the TIDEP-01012
%                       SSD (..._data.bin, ..._idx.bin)
%                    b) The json-Files describing the acquisition/settings
%                       with the ending .mmwave.json and .setup.json.
%
%
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang, ETH Zürich (02nd August 2022)

function [path2calibMatDir] = cascade_MIMO_00A_calibration_phaseMismatch(path2file, trg_rng, frame_ids)

dataPlatform = 'TDA2';

path2calibMatDir = fullfile(path2file,'00_Calibration_Files');
input_path_radar = fullfile(path2file,'01_Raw_Radar_Data');

if ~isfolder(path2calibMatDir)
    mkdir(path2calibMatDir);
end

%parameter file name
pathGenParaFile = fullfile(path2calibMatDir,'generateCalibrationMatrix_param.m');

%generate parameter file
parameter_file_gen_antennaCalib_json(input_path_radar, pathGenParaFile, dataPlatform);

for fi = 1:length(frame_ids)
    %calibrateValFileNameSave: file name to save calibration results. This file
    %will be saved in "dataFolder_calib" after running calibration
    calibrateValFileNameSave = fullfile(path2calibMatDir,sprintf('phaseMismatchCalibration_%07d.mat',frame_ids(fi)));

    %important to clear the same.m file, since Matlab does not clear cache
    %automatically
    clear(pathGenParaFile);

    % Create Calibration Matrix Object
    genCalibrationMatrixObj = genCalibrationMatrixCascade(...
                                'pfile', pathGenParaFile,...
                                'calibrateFileName',input_path_radar,...
                                'targetRange', trg_rng);

    [fileIdx_unique] = getUniqueFileIdx(input_path_radar);
    [fileNameStruct]= getBinFileNames_withIdx(input_path_radar, fileIdx_unique{1});    
    genCalibrationMatrixObj.binDataFile = fileNameStruct;

    if length(genCalibrationMatrixObj.TxToEnable)< 12
        %it is important to know that all 12 TX needs to be turned on in the MIMO mode to generate a correct calibration matrix. 
        error('This data set cannot be used for calibration, all 12 channels should be enabled');
    end

    genCalibrationMatrixObj.frameIdx = frame_ids(fi);

    calibResult = dataPath(genCalibrationMatrixObj);
    RangeMat = calibResult.RangeMat;
    targetRange_est = (floor(mean(RangeMat(:))/genCalibrationMatrixObj.calibrationInterp))...
        *genCalibrationMatrixObj.rangeResolution;

    fprintf('Target is estimated at range %.3fm (%cr=%0.3fm)\n',...
             targetRange_est,...
             char(916),...
             trg_rng-targetRange_est);

    if trg_rng-targetRange_est > 1
        fprintf('Last Valid Frame No. %d\n',frame_ids(fi)-1);
        fprintf('-------------------------------------------------\n');
        break
    end

    %just to make it compatible with old data
    params.Slope_MHzperus = genCalibrationMatrixObj.Slope_calib/1e12;
    params.Sampling_Rate_sps = genCalibrationMatrixObj.Sampling_Rate_sps;

    %save the calibration data
    save(calibrateValFileNameSave, 'calibResult','params');
end

end
