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

% parameter_gen_from_Jason.m
%
% function to generate paramsConfig from Json file
%input:
%   test_name: test data file name which contains the chirp configuration Json file
%   params: with some init value set before calling this function
%output:
%   params: output parameters
%generate paramsConfig from Json file

function params = parameter_gen_from_Jason(test_name, params)
%% choose which TX to used in TX TDM MIMO mode, all 16 RXs are enabled by default
% TI cascade board antenna configurations

platform = 'TI_4Chip_CASCADE';

%% fixed antenna ID and postion values for TI 4-chip cascade board. Should not be changed if user is based on TI board
% TI_Cascade_TX_position_azi = [0 4 8 12 16 20 24 28 32 9 10 11];%12 TX antenna azimuth position on TI 4-chip cascade EVM
% TI_Cascade_TX_position_ele = [0 0 0 0 0 0 0 0 0 1 4 6];%12 TX antenna elevation position on TI 4-chip cascade EVM
% %TI_Cascade_RX_position_azi = [0:7 39:42 50:53];%16 RX antenna azimuth position on TI 4-chip cascade EVM
% TI_Cascade_RX_position_azi = [0:3 11:14 50:53 46:49 ];
% TI_Cascade_RX_position_ele = zeros(1,16);%16 RX antenna elevation position on TI 4-chip cascade EVM
% TI_Cascade_RX_ID = [1 2 3 4 5 6 7 8 13 14 15 16 9 10 11 12]; %RX channel order on TI 4-chip cascade EVM

% New Format TDA
TI_Cascade_TX_position_azi = [11 10 9 32 28 24 20 16 12 8 4 0 ];%12 TX antenna azimuth position on TI 4-chip cascade EVM
TI_Cascade_TX_position_ele = [6 4 1 0 0 0 0 0 0 0 0 0];%12 TX antenna elevation position on TI 4-chip cascade EVM
TI_Cascade_RX_position_ele = zeros(1,16);%16 RX antenna elevation position on TI 4-chip cascade EVM
TI_Cascade_RX_position_azi = [ 11:14 50:53 46:49 0:3  ];
TI_Cascade_RX_ID = [13 14 15 16 1 2 3 4 9 10 11 12 5 6 7 8 ]; %RX channel order on TI 4-chip cascade EVM


TI_Cascade_Antenna_DesignFreq = 76.8; % antenna distance is designed for this frequency

speedOfLight = 3e8;

params.TI_Cascade_RX_ID = TI_Cascade_RX_ID;%RX channel IDs to use
D_RX = TI_Cascade_RX_position_azi(TI_Cascade_RX_ID); %RX azimuth antenna coordinate

params.D_RX = D_RX; %RX channel azimuth location in the unit of half wavelength
params.Rx_Elements_To_Capture = 1:16;  %All Rx
params.numRX = length(D_RX);

%% find *.json file in the same folder of the test data; if json files in the folder, this code will NOT work, please make sure only one exist
f=dir(fullfile(test_name,'*.mmwave.json'));
if length(f)~=1
    error('Unknown parameter file or too many chirpProfiles_*.json file!!');
else
    
    disp(['paramFile= ' test_name f.name]);
    paramFile = fullfile(test_name,f.name);
    
    %go to the test data folder to run the .json file and read the chirp parameters
    %associated with the data
    params_chirp = JsonParser(paramFile);
    
    
    
    %% chirp/profile parameters
    nchirp_loops = params_chirp.DevConfig(1).AdvFrame.SubFrame(1).NumChirpLoops;
    Num_Frames = params_chirp.DevConfig(1).AdvFrame.NumFrames;
    params.Start_Freq_GHz			=  params_chirp.DevConfig(1).Profile.StartFreq;								% For PG1 silicon, please avoid 76-77GHz band
    params.Slope_MHzperus			=  params_chirp.DevConfig(1).Profile.FreqSlope ;                              % MHz/us
    params.Idle_Time_us            =	params_chirp.DevConfig(1).Profile.IdleTime;                              % us
    params.Adc_Start_Time_us		=	params_chirp.DevConfig(1).Profile.AdcStartTime;                              % us
    params.Ramp_End_Time_us		=	params_chirp.DevConfig(1).Profile.RampEndTime;                             % us
    params.Sampling_Rate_ksps		=	params_chirp.DevConfig(1).Profile.SamplingRate;							% ksps
    params.Samples_per_Chirp		=	params_chirp.DevConfig(1).Profile.NumSamples;    						% Number of samples per chirp
    % Frame config
    params.nchirp_loops			=	nchirp_loops;							% Number of chirps per frame
    params.Num_Frames				=	Num_Frames;								% 42 frames gives 4GB
    params.Dutycycle               =   0.5;                            % (ON duration)/(ON+OFF duration)
    params.Chirp_Duration_us       =   (params.Ramp_End_Time_us+params.Idle_Time_us); % us
    params.NumberOfSamplesPerChannel = params.Samples_per_Chirp * nchirp_loops * params.NumAnglesToSweep *params.Num_Frames;
    
    
    centerFrequency = params.Start_Freq_GHz+(params.Samples_per_Chirp/params.Sampling_Rate_ksps*params.Slope_MHzperus)/2;
    d = 0.5*centerFrequency/TI_Cascade_Antenna_DesignFreq;
    params.d_BF = d;
    
    
    %advanced frame config
    params.numSubFrames = 1;
    %paramters for subframe 1
    if params.Chirp_Frame_BF == 0 % frame based
        params.SF1ChirpStartIdx = 1;
        params.SF1NumChirps = 1;
        params.SF1NumLoops = nchirp_loops;
        %multiply 200 to convert the value to be programed to the register %example:2000000=10ms
        params.SF1BurstPeriodicity = (params.Ramp_End_Time_us +params.Idle_Time_us)...
            *nchirp_loops /params.Dutycycle*200;              %example:2000000=10ms
        params.SF1ChirpStartIdxOffset = 1;
        params.SF1NumBurst = params.NumAnglesToSweep;
        params.SF1NumBurstLoops = 1;
        params.SF1SubFramePeriodicity = params.SF1BurstPeriodicity*params.NumAnglesToSweep;%50ms
    elseif params.Chirp_Frame_BF == 1 % chirp based
        params.SF1ChirpStartIdx = 1;
        params.SF1NumChirps = params.NumAnglesToSweep;
        params.SF1NumLoops = nchirp_loops;
        %multiply 200 to convert the value to be programed to the register  %example:2000000=10ms
        params.SF1BurstPeriodicity = (params.Ramp_End_Time_us +params.Idle_Time_us)...
            *nchirp_loops /params.Dutycycle*200 * params.NumAnglesToSweep ;
        params.SF1ChirpStartIdxOffset = 1;
        params.SF1NumBurst = 1;
        params.SF1NumBurstLoops = 1;
        params.SF1SubFramePeriodicity = params.SF1BurstPeriodicity;
        
    end
    
    params.Frame_Repetition_Period_ms = params.SF1SubFramePeriodicity/200/1000;
    
    params.calibrationInterp = 5;
    params.phaseCalibOnly = 1;
    
    %% algorithm parameters
    params.ApplyRangeDopplerWind = 1;
    params.rangeFFTSize = 2^ceil(log2(params.Samples_per_Chirp));
    %% derived parameters
    chirpRampTime       = params.Samples_per_Chirp/(params.Sampling_Rate_ksps/1e3);
    chirpBandwidth      = params.Slope_MHzperus(1) * chirpRampTime; % Hz
    rangeResolution     = speedOfLight/2/(chirpBandwidth*1e6);
    params.rangeBinSize        = rangeResolution*params.Samples_per_Chirp/params.rangeFFTSize;
end


end