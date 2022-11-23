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


%% Calculates phase shifter values for beam steering
% input: 
%      b_ang: steering angles
%      PS_actual: acutal measured phase value from calibration file,
%      different from programmed phase.
%      D_TX:  TX antenna array position in the unit of half lamda
%      d: = 0.5*centerFrequency/TI_Cascade_Antenna_DesignFreq 

% output: 
%      radar_data_frame_total_sub1: raw adc data


function [slope,PS_Tx,PS_forAoA]=BeamSteerPhaseCalc(b_ang,PS_actual, D_TX,d)

numTX = size(PS_actual,1);
PS_allowed=0:5.625:360-5.625;						%allowed phase shifts
actualPS10=-PS_actual;
actualPS10=[zeros(numTX,1) actualPS10];				%rearranged phase shifter LUT
PS_Tx=zeros(numTX,length(b_ang));						%Matrix to store actual phase shift values to be given to phase shifter
PS_forAoA=zeros(numTX,length(b_ang));					%Matrix to store actual phase shift that happens as maintained by LUT
slope=zeros(length(b_ang),numTX);						%Ideally calculated phase shifts

for t=1:length(b_ang)
    slope(t,:)= wrapTo360(sin(b_ang(t)*2*pi*d/180)*D_TX*180);
    for u=1:2:numTX		%finds nearest higher phase shift
        temp=actualPS10(u,:);
        diff=find((temp-slope(t,u))>=0);
        if(isempty(diff))
            diff=1;
        end
        PS_forAoA(u,t)=temp(diff(1));
        PS_Tx(u,t)=PS_allowed(diff(1));
    end
    for u=2:2:numTX	%finds nearest lower phase shift
        temp=actualPS10(u,:);
        diff=find((temp-slope(t,u))<=0);
        PS_forAoA(u,t)=temp(diff(end));
        PS_Tx(u,t)=PS_allowed(diff(end));
    end
end
end