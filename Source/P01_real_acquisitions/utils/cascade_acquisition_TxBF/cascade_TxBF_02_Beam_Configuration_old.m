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

%Main file for TX beamforming towards desired angles

function [ ] = cascade_TxBF_02_Beam_Configuration(phaseShiftCalibfile, phaseMismatchCalibfile, paramFile, rst_dll_path)

LOAD_CALIBRRATION  = 1;

%phase Shifter calibration file (TxBF)
%phaseShiftCalibfile = 'Z:\Publications\01_Paper\2022\xx_JAG_MIMO-SAR_ABOetal\Matlab_Functions\D00_sample_data\real\TxBF_A77_phaseShiftCalibration_LR230\00_Calibration_Files\phaseShifterCalibration.mat';
%load phase Mismatch calibration vector for TX channels (MIMO)
%phaseMismatchCalibfile = 'Z:\Publications\01_Paper\2022\xx_JAG_MIMO-SAR_ABOetal\Matlab_Functions\D00_sample_data\real\MIMO_A77_phaseMismatchCalibration_LR230\00_Calibration_Files\phaseMismatchCalibration.mat';
%chirp config parameters to use
%paramFile = 'chirpProfile_TxBF_LR230';

cmdStr = sprintf('paramsConfig = %s();',paramFile);
eval(cmdStr);

% Initialize Radarstudio .NET connection
if ~exist('rst_dll_path','var')
    rst_dll_path = 'C:\ti\mmwave_studio_03_00_00_14\mmWaveStudio\Clients\RtttNetClientController\RtttNetClientAPI.dll';
end

ErrStatus = Init_RSTD_Connection(rst_dll_path);
if (ErrStatus ~= 30000)
    disp('Error inside Init_RSTD_Connection');
    return;
end

%loads phase shifter calibrations. Measured zero angle data in DEGREES,
load(phaseShiftCalibfile);
PS_actual = (squeeze((Ph(:,1,:))))';
PS_actual(PS_actual > 0)  =  PS_actual(PS_actual>0) - 360;

if (LOAD_CALIBRRATION == 1)
    
    load(phaseMismatchCalibfile);
    Rxcal = calibResult.RxMismatch;
    CorrectPhase = calibResult.TxMismatch;
else
    CorrectPhase=zeros(3,4);
end


% if (numel(Rxcal) ~= length(paramsConfig.Rx_Elements_To_Capture))
%     disp('Number of Rx elements in cal file not matching Rx MIMO array');
%     return
% end
clear AngleMat
clear RangeMat
clear ValueMat

%calculate the phase values for each angle to be steered
[slope,PS_Tx,PS_forAoA]=BeamSteerPhaseCalc(paramsConfig.anglesToSteer,PS_actual(paramsConfig.Tx_Ant_Arr_BF,:),paramsConfig.D_TX_BF, paramsConfig.d_BF);
for iter=1:length(paramsConfig.anglesToSteer)
    
    TXphaseMat_test = zeros(1,12); % 12 is the number of TX channels
    TXphaseMat_test(paramsConfig.Tx_Ant_Arr_BF) = PS_Tx(:,iter);
    TXphaseMat = reshape(TXphaseMat_test, [3 4]); % 3 the numTX channel per device; 4 the number of devices
    %disp(['Beam Angle from normal ' num2str(paramsConfig.anglesToSteer)]);
    
    % add the initial calibration phase value
    Tx1_Ph_Deg(iter,:)=wrapTo360(CorrectPhase(1,:)+TXphaseMat(1,:));
    Tx2_Ph_Deg(iter,:)=wrapTo360(CorrectPhase(2,:)+TXphaseMat(2,:));
    Tx3_Ph_Deg(iter,:)=wrapTo360(CorrectPhase(3,:)+TXphaseMat(3,:));
end


%config each device with cprresponding chirp parameters and  phase values
for ar_device = 1:4
    ar_device
    if(paramsConfig.RadarDevice(ar_device)==1)
        
        [ErrStatus] ...
            = Sensor_AdvanceFrame_BF(ar_device,paramsConfig,Tx1_Ph_Deg,Tx2_Ph_Deg,Tx3_Ph_Deg);
        
        if (ErrStatus ~= 30000)
            disp('Error inside Sensor_Config_Function');
            return;
        end
    end
end

end