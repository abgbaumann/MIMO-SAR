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

%max range: 300 meters

function [params] = chirpProfile_TxBF_LR230_30deg_1chirp()



% TI cascade board antenna configurations
platform = 'TI_4Chip_CASCADE';

%% fixed antenna ID and postion values for TI 4-chip cascade board. Should not be changed if user is based on TI board
TI_Cascade_TX_position_azi = [11 10 9 32 28 24 20 16 12 8 4 0 ];%12 TX antenna azimuth position on TI 4-chip cascade EVM
TI_Cascade_TX_position_ele = [6 4 1 0 0 0 0 0 0 0 0 0];%12 TX antenna elevation position on TI 4-chip cascade EVM
TI_Cascade_Antenna_DesignFreq = 76.8; % antenna distance is designed for this frequency
speedOfLight = physconst('LightSpeed');
params.Tx_Ant_Arr_BF = [ 12:-1:4]; %TX channel IDs to use for beamforming; all 16 RX channels are enabled by default
params.D_TX_BF = TI_Cascade_TX_position_azi(params.Tx_Ant_Arr_BF); %TX azimuth antenna coordinates in unit of half lamda

% params.D_RX = D_RX; %RX channel azimuth location in the unit of half wavelength
params.RadarDevice=[1 1 1 1];			%set 1 all the time
params.Rx_Elements_To_Capture = 1:16;   %All Rx


%% define what angles to steer in TX beamforming mode
params.anglesToSteer=[-30:2:30];   % angles to steer in TX beamforming mode, in unit of degrees
params.NumAnglesToSweep = length(params.anglesToSteer);


%% chirp/profile parameters
nchirp_loops                        =   1;
Num_Frames                          =   100;
f_acq                               =   10;              % acquistion frequency [Hz] = [1/s]
t_acq_rep                           =   1 / f_acq;      % frame itnerval [s]
params.Start_Freq_GHz			    =	77;				% starting fequency for chirp, make sure the entire BW is within 76~77 or 77~81
params.Slope_MHzperus			    =   4.9728;         % MHz/us
params.Idle_Time_us                 =	5;              % us
params.Tx_Start_Time_us             =   0;              % us
params.Adc_Start_Time_us		    =	10;             % us
params.Ramp_End_Time_us		        =	80;             % us
params.Sampling_Rate_ksps		    =	7650;			% ksps
params.Samples_per_Chirp		    =	512;    	    % Number of samples per chirp
params.Rx_Gain_dB				    =	30;			    % dB
% Frame config
params.nchirp_loops			        =	nchirp_loops;   % Number of chirps per frame
params.Num_Frames				    =	Num_Frames;		% number of frames to collect data
params.Chirp_Duration_us            =   (params.Ramp_End_Time_us+params.Idle_Time_us); % us
SingleAngleDuration_s               =   params.Chirp_Duration_us * 10^-6 * nchirp_loops; % [s]
SingleFrameDuration_s               =   SingleAngleDuration_s * params.NumAnglesToSweep; % [s]

if SingleFrameDuration_s <= t_acq_rep
    params.Dutycycle                    =   SingleFrameDuration_s / t_acq_rep;
else
    scale = 0.999;
    warning(sprintf('Acquisition Frequency is too high. Set to maximum available frequency. (1/%d Hz)',round(1/SingleFrameDuration_s*scale,0)))
    params.Dutycycle                    =   SingleFrameDuration_s / SingleFrameDuration_s * scale;
end

%params.Dutycycle                    =   SingleAngleDuration_s / SingleFrameDuration_s;
%params.Dutycycle                    =   0.5;            % (ON duration)/(ON+OFF duration)
params.NumberOfSamplesPerChannel    =   params.Samples_per_Chirp * nchirp_loops * params.NumAnglesToSweep *params.Num_Frames; %number of ADC samples received per channel. this value is used in HSDC for data capture


%d = 0.5*actual wavelength/wavelength for antenna design = 0.5 * actual center frequency/board antenna design frequency
centerFrequency = params.Start_Freq_GHz+(params.Samples_per_Chirp/params.Sampling_Rate_ksps*params.Slope_MHzperus)/2;
d = 0.5*centerFrequency/TI_Cascade_Antenna_DesignFreq;
params.d_BF = d;


%advanced frame config
params.Chirp_Frame_BF = 0; % 1 - chirp based beam steering, 0 - frame based beam steering
params.numSubFrames = 1;
%paramters for subframe 1
if params.Chirp_Frame_BF == 0 % frame based
    params.SF1ChirpStartIdx = 0; %SF1 Start index of the first chirp in this sub frame
    params.SF1NumChirps = 1; % SF1 Number Of unique Chirps per burst
    params.SF1NumLoops = nchirp_loops;% SF1 Number Of times to loop through the unique chirps in each burst
    %multiply 200 to convert the value to be programed to the register %example:2000000=10ms
    params.SF1BurstPeriodicity = (params.Ramp_End_Time_us + params.Idle_Time_us)...
        * nchirp_loops / params.Dutycycle * 200;              %example:2000000=10ms
    params.SF1ChirpStartIdxOffset = 1; % SF1 Chirps Start Idex Offset
    params.SF1NumBurst = params.NumAnglesToSweep; % SF1 Number Of Bursts constituting this sub frame
    params.SF1NumBurstLoops = 1; % SF1 Number Of Burst Loops
    params.SF1SubFramePeriodicity = params.SF1BurstPeriodicity*params.NumAnglesToSweep;%SF1 SubFrame Periodicity
elseif params.Chirp_Frame_BF == 1 % chirp based
    params.SF1ChirpStartIdx = 0;
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

%frame repitition time in unit of ms. THis based on params.SF1SubFramePeriodicity
params.Frame_Repetition_Period_ms = params.SF1SubFramePeriodicity/200/1000;


%% algorithm parameters
params.ApplyRangeDopplerWind = 1;
params.rangeFFTSize = 2^ceil(log2(params.Samples_per_Chirp));


%% derived parameters
chirpRampTime       = params.Samples_per_Chirp/(params.Sampling_Rate_ksps/1e3);
chirpBandwidth      = params.Slope_MHzperus(1) * chirpRampTime; % Hz
rangeResolution     = speedOfLight/2/(chirpBandwidth*1e6);
params.rangeBinSize = rangeResolution*params.Samples_per_Chirp/params.rangeFFTSize;


end



