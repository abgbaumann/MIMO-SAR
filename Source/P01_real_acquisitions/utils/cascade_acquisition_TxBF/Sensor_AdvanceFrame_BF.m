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


%% beamforming using advanced frame configuration, support both chirp based and frame based beam steering
% input: 
%      ar_device: device ID
%      params: configuration parameters
%      Tx1_Ph_Deg:  phase value for TX1 channels for all 4 devices for all desired angles
%      Tx2_Ph_Deg:  phase value for TX2 channels for all 4 devices for all desired angles
%      Tx3_Ph_Deg:  phase value for TX3 channels for all 4 devices for all desired angles


% output: 
%      ErrStatus: error status


function  [ErrStatus ] = Sensor_AdvanceFrame_BF(ar_device,params,Tx1_Ph_Deg,Tx2_Ph_Deg,Tx3_Ph_Deg)

%convert to angle step size of 5.625
Tx1_Ph_Deg = floor(Tx1_Ph_Deg / 5.625);
Tx2_Ph_Deg = floor(Tx2_Ph_Deg / 5.625);
Tx3_Ph_Deg = floor(Tx3_Ph_Deg / 5.625);

%number of TXs for TDM MIMO
Tx_Ant_Arr_BF = params.Tx_Ant_Arr_BF;

NumAnglesToSweep = params.NumAnglesToSweep;
Num_Frames = params.Num_Frames;



%Misc parameters
% wait_time				=	1;                              %In Seconds
RadarDeviceId			=	[1, 2, 4, 8];
TriggerSelect_Arr		=	[1, 2, 2, 2];                   %1: Software trigger, 2: Hardware trigger

%%%%%%%%%%%%%-Sensor Configuration%%%%%%%%%%%%-
%Profile Configuration
%Syntax:
% Int32 ar1.ProfileConfig_mult(UInt16 RadarDeviceId, UInt16 profileId, Single startFreqConst, Single idleTimeConst, Single adcStartTimeConst, Single rampEndTime, UInt32 tx1OutPowerBackoffCode, UInt32 tx2OutPowerBackoffCode, UInt32 tx3OutPowerBackoffCode, UInt16 tx1PhaseShifter, UInt16 tx2PhaseShifter, UInt16 tx3PhaseShifter, Single freqSlopeConst, Single txStartTime, UInt16 numAdcSamples, UInt16 digOutSampleRate, Char hpfCornerFreq1, Char hpfCornerFreq2, Char rxGain) - Profile configuration API which defines chirp profile parameters
% _I_ UInt16	RadarDeviceId	 - Radar Device Id
% _I_ UInt16	profileId	 - Chirp Profile Id [0 to 3]
% _I_ Single	startFreqConst	 - Chirp Start Frequency in GHz
% _I_ Single	idleTimeConst	 - Chirp Idle Time in 탎
% _I_ Single	adcStartTimeConst	 - Chirp ADC Start Time in 탎
% _I_ Single	rampEndTime	 - Chirp Ramp End Time in 탎
% _I_ UInt32	tx1OutPowerBackoffCode	 - TX1 channel Power Backoff in dB
% _I_ UInt32	tx2OutPowerBackoffCode	 - TX2 channel Power Backoff in dB
% _I_ UInt32	tx3OutPowerBackoffCode	 - TX3 channel Power Backoff in dB
% _I_ UInt16	tx1PhaseShifter	 - TX1 channel Phase Shifter Value in deg
% _I_ UInt16	tx2PhaseShifter	 - TX2 channel Phase Shifter in deg
% _I_ UInt16	tx3PhaseShifter	 - TX3 channel Phase Shifter in deg
% _I_ Single	freqSlopeConst	 - Chirp Frequency Slope in MHz/탎
% _I_ Single	txStartTime	 - TX Start Time in 탎
% _I_ UInt16	numAdcSamples	 - RX Number of Adc Samples
% _I_ UInt16	digOutSampleRate	 - RX Sampling Rate in ksps
% _I_ Char	hpfCornerFreq1	 - RX HPF1 corner frequency, 0x00-175 kHz, 0x01-235 kHz, 0x02-350 kHz, 0x03-700 kHz
% _I_ Char	hpfCornerFreq2	 - RX HPF2 corner frequency, 0x00-350 kHz, 0x01-700 kHz, 0x02-1.4 MHz, 0x03-2.8 MHz
% _I_ Char	rxGain	 - RX Gain in dB

%set up one profile for TX BF
%profile for TX beamforming
Lua_String = [];
profileID = 0;

Lua_String = sprintf('ar1.ProfileConfig_mult(%d, %d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %d, %d, 2, 1, %d)', ...
    RadarDeviceId(ar_device), profileID, params.Start_Freq_GHz, params.Idle_Time_us, params.Adc_Start_Time_us, ...
    params.Ramp_End_Time_us, 0, 0, 0, 0, 0,0, ...
    params.Slope_MHzperus, params.Tx_Start_Time_us, params.Samples_per_Chirp,...
    params.Sampling_Rate_ksps, params.Rx_Gain_dB);
ErrStatus =RtttNetClientAPI.RtttNetClient.SendCommand(Lua_String);
disp('HPF cutoff increased to 350 and 700');
if (ErrStatus ==30000)
    disp('Profile Configuration successful for BF')
else
    disp('Profile Configuration failed for BF')
    return;
end


%Chirp Configuration
        %Syntax:
        % Int32 ar1.ChirpConfig_mult(UInt16 RadarDeviceId, UInt16 chirpStartIdx, UInt16 chirpEndIdx, UInt16 profileId, Single startFreqVar, Single freqSlopeVar, Single idleTimeVar, Single adcStartTimeVar, UInt16 tx1Enable, UInt16 tx2Enable, UInt16 tx3Enable) - Chirp configuration API which defines which profile is to be used for each chirp in a frame
        % _I_ UInt16	RadarDeviceId	 - Radar Device Id
        % _I_ UInt16	chirpStartIdx	 - First Chirp Start Index number
        % _I_ UInt16	chirpEndIdx	 - Last chirp Index number
        % _I_ UInt16	profileId	 - Chirp Configured profileId
        % _I_ Single	startFreqVar	 - Chirp start frequency var in MHz
        % _I_ Single	freqSlopeVar	 - frequency Slope Var in MHz/탎
        % _I_ Single	idleTimeVar	 - Idle Time Var in 탎
        % _I_ Single	adcStartTimeVar	 - ADC Start Time Var in 탎
        % _I_ UInt16	tx1Enable	 - tx1 channel
        % _I_ UInt16	tx2Enable	 - tx2 channel
        % _I_ UInt16	tx3Enable	 - tx3 channel
%set up chirp configuration for TXBF
cnt = 1;
profileID = 0;
for ChirpCfg =(0:NumAnglesToSweep-1)
    TXenableMat=zeros(3,4);
    TXenableMat(Tx_Ant_Arr_BF) = 1;
    Tx1_Enable = TXenableMat(1,:);
    Tx2_Enable = TXenableMat(2,:);
    Tx3_Enable = TXenableMat(3,:);
    
    Lua_String  = sprintf('ar1.ChirpConfig_mult(%d, %d, %d, %d, 0, 0, 0, 0, %d, %d, %d)', RadarDeviceId(ar_device), ...
        ChirpCfg, ChirpCfg,profileID, Tx1_Enable(ar_device), Tx2_Enable(ar_device), Tx3_Enable(ar_device));
    ErrStatus   =   RtttNetClientAPI.RtttNetClient.SendCommand(Lua_String);
%     if (ErrStatus ==30000)
%         disp('Chirp config successful')
%     else
%         disp('Chirp config failed')
%         return;
%     end
    Lua_String = [];
    Lua_String  = sprintf('ar1.SetPerChirpPhaseShifterConfig_mult(%d, %d, %d, %d, %d, %d)', RadarDeviceId(ar_device),...
        ChirpCfg, ChirpCfg, Tx1_Ph_Deg(cnt, ar_device), Tx2_Ph_Deg(cnt,ar_device), Tx3_Ph_Deg(cnt,ar_device));
    ErrStatus   =   RtttNetClientAPI.RtttNetClient.SendCommand(Lua_String);
%     if (ErrStatus ==30000)
%         disp('Phase shifter per Chirp config successful')
%     else
%         disp('Phase shifter per Chirp config failed')
%         return;
%     end
    cnt = cnt + 1;
    
    %pause(wait_time);
end

if (ErrStatus ==30000)
    disp('Chirp and Phase shifter per Chirp config successful')
else
    disp('Chirp or Phase shifter per Chirp config failed')
    return;
end

%advanced frame config
numSubFrames = params.numSubFrames;
%paramters for subframe 1
SF1ChirpStartIdx = params.SF1ChirpStartIdx;
SF1NumChirps = params.SF1NumChirps;
SF1NumLoops = params.SF1NumLoops;
SF1BurstPeriodicity = params.SF1BurstPeriodicity; 
SF1ChirpStartIdxOffset = params.SF1ChirpStartIdxOffset;
SF1NumBurst = params.SF1NumBurst;
SF1NumBurstLoops = params.SF1NumBurstLoops;
SF1SubFramePeriodicity = params.SF1SubFramePeriodicity;

Lua_String = [];
Lua_String = sprintf('ar1.AdvanceFrameConfig_mult(%d, %d, 1536, 0, %d, %d, %d, %d, %d, %d, %d, %d, 0, %d, %d, %d, %d, %d,%d, %d, %d, 0, 0, 1, 128, 8000000, 0, 1, 1, 8000000, 0, 0, 1, 128,8000000, 0, 1, 1, 8000000, %d, %d, 0, 0, 128, 256, 1, 128, 256, 1, 128,1, 1, 128, 1, 1)',...
    RadarDeviceId(ar_device), numSubFrames, SF1ChirpStartIdx, SF1NumChirps, SF1NumLoops, SF1BurstPeriodicity, ...
        SF1ChirpStartIdxOffset, SF1NumBurst, SF1NumBurstLoops, SF1SubFramePeriodicity, ...
    SF1ChirpStartIdx, SF1NumChirps, SF1NumLoops, SF1BurstPeriodicity,SF1ChirpStartIdxOffset,SF1NumBurst, ...
    SF1NumBurstLoops ,SF1SubFramePeriodicity, Num_Frames, TriggerSelect_Arr(ar_device));

ErrStatus =RtttNetClientAPI.RtttNetClient.SendCommand(Lua_String);
if (ErrStatus==30000)  %n*2 for the work around to handle an issue in the TSW14J56 FW
    disp('Frame Config successful');
else
    disp('Frame Conig failed');
    return;
end
%pause(wait_time);
disp(['AR Device ' num2str(ar_device) ' Configuration Successful']);


return;
end