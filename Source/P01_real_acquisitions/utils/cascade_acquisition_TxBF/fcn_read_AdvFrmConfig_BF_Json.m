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


%% read raw adc data TX beamforming from advaced configuration

% input:
%      BF_data_dest: data file path to be read and analyzed
%      params:  configuration parameters
% output:
%      radar_data_frame_total_sub1: raw adc data

function [adcData1Complex] = fcn_read_AdvFrmConfig_BF_Json(fileNameCascade, params)

numChirpPerLoop = params.nchirp_loops;
N_angles = params.NumAnglesToSweep;
numSamplePerChirp = params.Samples_per_Chirp;
numRX = params.numRX;
numDevices = 4;
frameIdx =  params.frameId;
numRXPerDevice = numRX/numDevices;

dataFolder =fileNameCascade.dataFolderName;
fileFullPath_master = fullfile(dataFolder,fileNameCascade.master);
fileFullPath_slave1 = fullfile(dataFolder,fileNameCascade.slave1);
fileFullPath_slave2 = fullfile(dataFolder,fileNameCascade.slave2);
fileFullPath_slave3 = fullfile(dataFolder,fileNameCascade.slave3);


if (params.Chirp_Frame_BF)
    numLoops = numChirpPerLoop;
    numChirpPerLoop = N_angles;
    [radar_data_Rxchain_master] = readBinFile(fileFullPath_master, params.Chirp_Frame_BF, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices);
    [radar_data_Rxchain_slave1] = readBinFile(fileFullPath_slave1, params.Chirp_Frame_BF, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices);
    [radar_data_Rxchain_slave2] = readBinFile(fileFullPath_slave2, params.Chirp_Frame_BF, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices);
    [radar_data_Rxchain_slave3] = readBinFile(fileFullPath_slave3, params.Chirp_Frame_BF, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices);
    adcData1Complex(:,:,1:4,:) = permute(radar_data_Rxchain_master, [2 4 1 3]);
    adcData1Complex(:,:,5:8,:) = permute(radar_data_Rxchain_slave1, [2 4 1 3]);
    adcData1Complex(:,:,9:12,:) = permute(radar_data_Rxchain_slave2, [2 4 1 3]);
    adcData1Complex(:,:,13:16,:) = permute(radar_data_Rxchain_slave3, [2 4 1 3]);
    
else
    numLoops = N_angles;
    [radar_data_Rxchain_master] = readBinFile(fileFullPath_master, params.Chirp_Frame_BF, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices);
    [radar_data_Rxchain_slave1] = readBinFile(fileFullPath_slave1, params.Chirp_Frame_BF, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices);
    [radar_data_Rxchain_slave2] = readBinFile(fileFullPath_slave2, params.Chirp_Frame_BF, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices);
    [radar_data_Rxchain_slave3] = readBinFile(fileFullPath_slave3, params.Chirp_Frame_BF, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices);
    adcData1Complex(:,:,1:4,:) = permute(radar_data_Rxchain_master, [2 3 1 4]);
    adcData1Complex(:,:,5:8,:) = permute(radar_data_Rxchain_slave1, [2 3 1 4]);
    adcData1Complex(:,:,9:12,:) = permute(radar_data_Rxchain_slave2, [2 3 1 4]);
    adcData1Complex(:,:,13:16,:) = permute(radar_data_Rxchain_slave3, [2 3 1 4]);
    
end



end



function [adcData1Complex] = readBinFile(fileFullPath, Chirp_Frame_BF, frameIdx,numSamplePerChirp,numChirpPerLoop,numLoops, numRXPerDevice, numDevices)
Expected_Num_SamplesPerFrame = numSamplePerChirp*numChirpPerLoop*numLoops*numRXPerDevice*2;
fp = fopen(fileFullPath, 'r');
fseek(fp,(frameIdx-1)*Expected_Num_SamplesPerFrame*2, 'bof');
adcData1 = fread(fp,Expected_Num_SamplesPerFrame,'uint16');
neg             = logical(bitget(adcData1, 16));
adcData1(neg)    = adcData1(neg) - 2^16;
%%
adcData1 = adcData1(1:2:end) + sqrt(-1)*adcData1(2:2:end);
adcData1Complex = reshape(adcData1, numRXPerDevice, numSamplePerChirp, numChirpPerLoop, numLoops);



fclose(fp);
end