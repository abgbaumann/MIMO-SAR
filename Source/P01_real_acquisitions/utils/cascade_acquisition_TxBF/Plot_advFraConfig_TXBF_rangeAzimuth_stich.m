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


% show range/azimuth heat map within scanned field of by stiching multiple
% scans together using RX array

% input:
%      BF_data: raw ADC data to be analyzed
%      Rx_Ant_Arr: RX channel ID order to be read in sequence
%      params:  configuration parameters
%      SweepAngles: angles that swept/steered in this data
%      BF_MIMO_ref: phase calibration vector for RX channels
%      startRangeInd: start range index to be included
%      lastRangeIndToThrow: number of last range bins to be excluded in the plot
% output:
%      range_angle_stich: stiched range azimuth heat map

function [ range_angle_stich, y_axis, x_axis]=...
    Plot_advFraConfig_TXBF_rangeAzimuth_stich(BF_data,Rx_Ant_Arr,params,SweepAngles,BF_MIMO_ref, startRangeInd, lastRangeIndToThrow)

D_RX = params.D_RX;
nfft3=256;											% 3rd FFT interpolation
Samples_per_Chirp = params.Samples_per_Chirp;
rangeFFTSize = 2^ceil(log2(Samples_per_Chirp));
d = params.d_BF;


Sampling_Rate_ksps = params.Sampling_Rate_ksps;
Slope_MHzperus = params.Slope_MHzperus;
num_sweep=length(SweepAngles);
DopplerFFTSize = 2^ceil(log2(size(BF_data,2)));
rangeResolution = 3e8/2/(Samples_per_Chirp/Sampling_Rate_ksps*Slope_MHzperus*1e9);
indices_1D = (startRangeInd: params.Samples_per_Chirp-lastRangeIndToThrow) ;
resolution = physconst('LightSpeed')/2/(params.Samples_per_Chirp/params.Sampling_Rate_ksps*1e3*Slope_MHzperus*1e6);


%window for range FFT
window_1D = repmat(hanning(Samples_per_Chirp),1, size(BF_data,2));
window_2D = repmat(hanning(DopplerFFTSize),1, size(BF_data,1));
window_2D = window_2D';

for Rxnum=1:length(Rx_Ant_Arr)
    for sweep=1:num_sweep
        radar_data_frame = squeeze(BF_data(:,:,Rxnum,sweep));        
       
        rangeFFT1(:,:,Rxnum,sweep) = fft(radar_data_frame.*window_1D,rangeFFTSize,1);
        rangeDopFFT1(:,:,Rxnum,sweep) = fftshift(fft(rangeFFT1(:,:,Rxnum,sweep).*window_2D,DopplerFFTSize,2),2);
        rangeDopFFT1(:,:,Rxnum,sweep) = rangeDopFFT1(:,:,Rxnum,sweep).*(exp(-1i*BF_MIMO_ref(Rxnum)*pi/180));
    end
end

%perform beamsteering towards the angle TX steering angles
rangeDopFFT1_zeroDop = squeeze(rangeDopFFT1(:,floor(DopplerFFTSize/2+1),:,:));
range_angle_stich = [];
for iangle = 1:num_sweep
    angleTX = SweepAngles(iangle);
    wx = sind(angleTX);
    a1_az = exp(1*j*2*pi*d*(D_RX*wx));
    
    for irange = 1:size(rangeDopFFT1_zeroDop,1)
        RX_data = squeeze(rangeDopFFT1_zeroDop(irange,:,iangle));
        range_angle_stich(irange, iangle) = a1_az*(RX_data'*RX_data)*a1_az';
    end
end


rangeDopFFT1_ave = mean(mean(abs(rangeDopFFT1),3),4);
range_scale = (0:size(rangeDopFFT1_ave,1)-1)* rangeResolution;

if length(SweepAngles) > 1
    sine_theta = sind(SweepAngles);
    cos_theta = sqrt(1-sine_theta.^2);
    range_1D = (1:params.Samples_per_Chirp);
    [R_mat, sine_theta_mat] = meshgrid(range_1D*resolution,sine_theta);
    [~, cos_theta_mat] = meshgrid(range_1D,cos_theta);
    x_axis = R_mat.*cos_theta_mat;
    y_axis = R_mat.*sine_theta_mat;
end

range_angle_stich = range_angle_stich';


end
