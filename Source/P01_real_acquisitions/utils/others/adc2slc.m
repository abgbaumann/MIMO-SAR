function [slcStatic,...
          slcDynamic,...
          x_along_axis,...
          y_cross_axis] = adc2slc(adcData, rngRes, filt_type, numAngles, ant_az_id)
% 
% Function to convert raw Analog-to-Digital-Converted data to Single Look
% Complex data and calculating the position in along range and cross range.
%
% -------------------------------------------------------------------------
% Input Parameters:
% 
% adcData:       Raw Analog-to-Digital-Converted Data (complex)
% rngRes:        Range Resolution
% filt_type:     Window Filter (i.e. 'hanning', 'blackman')
% numAngles:     Number of Angles, for SAR: numAngles = N_VA; (e.g. 86)
% ant_az_id:     Indices of Azimuth Antenna         
%
% -------------------------------------------------------------------------
% Output Parameters:
%
% slcStatic:     Static Single Look Complex
% slcDynamic:    Dynamic Single Look Complex (Doppler)
% x_along_axis:  Along Range Coordinates         
% y_cross_axis:  Cross Range Coordinates    
%
% -------------------------------------------------------------------------
% How to use:
%
% numRange = 512;               % Number of Range Bins
% numChirp = 32;                % Number of Chirps per Loop
% numRXA = 16;                  % Number of RX antennas
% numTXA = 12;                  % Number of TX antennas
% data_sz = [numRange,numChirp,numRXA,numTXA]; % Size of adcData
% adcData = complex(rand(data_sz),rand(data_sz)); % raw ADC data
% rngRes = 0.1;                 % Range Resolution
% filt_type = 'hanning';        % Filter Type
% ant_az_id = 1:numRXA*numTXA;  % Azimuth antenna IDs in Consecutive Order
% [slcS,slcD,x,y] = adc2slc(adcData, rngRes, filt_type, ant_az_id);
%
% -------------------------------------------------------------------------
% Copyright (C) 2021 Andreas Baumann-Ouyang, ETH Zürich

ploting_fig = 0;
[numRange, numChirp, numRXA, numTXA] = size(adcData);
%numAngles = 256; % for SAR resolution: numAngles = N_VA; (e.g. 86)
RangeBinDiscard_atEnd = 0; % Remove the last N range bins.
RangeBinDiscard_atStart = 0; % Remove the first N range bins.
filt_az = 0; % Apply Window Filtering on Azimuth

if ~exist('ant_az_id','var')
   ant_az_id = 1:numRXA*numTXA; 
end

if strcmp(filt_type,'hanning')
    filt = hanning(numRange);
elseif strcmp(filt_type,'blackman')
    filt = blackman(numRange);
elseif strcmp(filt_type,'taylor')
    % W.-D. Wirth. Radar Techniques Using Array Antennas. 2nd ed. 
    % The Institution of Engineering and Technology, 2013. 
    % ISBN: 978-1-84919-698-7.
    filt = taylorwin(numRange);
end

%% Range Processing
rngData = complex(zeros(numRange, numChirp, numRXA, numTXA));
for txi = 1:numTXA
    for rxi = 1:numRXA
        data = adcData(:,:,rxi,txi);
        data = bsxfun(@minus, data, mean(data));
        data = bsxfun(@times, data, filt);
        rngData(:,:,rxi,txi) = fft(data, numRange, 1);
    end
end

%% Doppler Processing
dopData = complex(zeros(numRange, numChirp, numRXA, numTXA));
for txi = 1:numTXA
    for rxi = 1:numRXA
        data = rngData(:,:,rxi,txi);
        data = bsxfun(@times, data, ones([1,numChirp]));
        dopData(:,:,rxi,txi) = fftshift(fft(data,numChirp,2),2);
    end
end

dopData = reshape(dopData, numRange, numChirp, numRXA*numTXA);

if numChirp>1
    for dpi = 1:numChirp
        deltaPhi = 2 * pi * ((dpi-1)- numChirp/2) / ...
                            (numTXA * numChirp);
        sig = dopData(:,dpi,:);
        for txi = 1:numTXA
            RXID = (txi-1)*numRXA+1:txi*numRXA;
            corVec = repmat(exp(-1i*(txi-1)*deltaPhi), ...
                            numRange, numRXA);
            dopData(:,dpi,RXID) = sig(:,RXID) .* corVec;
        end
    end
end

dopDataAz = dopData(:,:,ant_az_id);

%% Azimuth Processing
if filt_az
    if strcmp(filt_type,'hanning')
        filt = hanning(length(ant_az_id));
    elseif strcmp(filt_type,'blackman')
        filt = blackman(length(ant_az_id));
    elseif strcmp(filt_type,'taylor')
        % W.-D. Wirth. Radar Techniques Using Array Antennas. 2nd ed. 
        % The Institution of Engineering and Technology, 2013. 
        % ISBN: 978-1-84919-698-7.
        filt = taylorwin(length(ant_az_id));
    end
    dopDataAz = bsxfun(@minus, dopDataAz, mean(dopDataAz,3));
    for chirp_i = 1:numChirp
        dopDataAz(:,chirp_i,:) = bsxfun(@times, squeeze(dopDataAz(:,chirp_i,:)), filt');
    end
end
azData = fft(dopDataAz,numAngles,3);

%% Dynamic SLC
if numChirp>1
    ratio = 0.5;
    DopplerPower = sum(mean((abs(azData(:,:,:))),3),1);
    DopplerPower_noDC = DopplerPower([1: numChirp/2-1 numChirp/2+3:end]);
    [peakVal, ~] = max(DopplerPower_noDC);
    threshold = peakVal*ratio;
    indSel = find(DopplerPower_noDC >threshold);
    for ii = 1:length(indSel)
        if indSel(ii) > numChirp/2-1
            indSel(ii) = indSel(ii) + 3;
        end
    end
    slcDynamic = squeeze(sum(azData(:,indSel,:),2));
    slcDynamic = fftshift(slcDynamic,2);
    slcDynamic = flipud(slcDynamic');
else
    slcDynamic = complex(zeros(size(squeeze(azData(:,1,:)))))';
end

%% Static SLC
if numChirp>1
    idx = numChirp/2+1;
else
    idx = 1;
end
slcStatic = squeeze(sum(azData(:,idx,:),2));
slcStatic = fftshift(slcStatic,2);
slcStatic = flipud(slcStatic');

%% Coordinates
theta = asin(-2*((-numAngles/2:numAngles/2)/numAngles));
theta_c = theta(1:end-1) + diff(theta)/2;
range_c = double(1:numRange) * rngRes - rngRes/2;

[range_mat, theta_mat] = meshgrid(range_c,theta_c);
x_along_axis = range_mat.*cos(theta_mat);
y_cross_axis = range_mat.*sin(theta_mat);

%% Discard last Range Bins
N = size(x_along_axis,2);
idx = RangeBinDiscard_atStart+1:N-(RangeBinDiscard_atEnd+1);
x_along_axis = x_along_axis(:,idx);
y_cross_axis = y_cross_axis(:,idx);
slcDynamic = slcDynamic(:,idx);
slcStatic = slcStatic(:,idx);

%% Ploting
if ploting_fig
    if numChirp>1
        amp = 20*log10(abs(slcDynamic));
        pha = angle(slcDynamic);
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(2,1,1);
        c=surf(y_cross_axis,x_along_axis,amp);
        view(0,90);
        axis equal
        set(c,'EdgeColor','none');
        title('Amplitude (Dynamic)');
        xlim([min(y_cross_axis(:)),max(y_cross_axis(:))]);
        ylim([min(x_along_axis(:)),max(x_along_axis(:))]);
        colorbar;
        subplot(2,1,2);
        c=surf(y_cross_axis,x_along_axis,pha);
        view(0,90);
        axis equal
        set(c,'EdgeColor','none');
        title('Phase (Dynamic)');
        xlim([min(y_cross_axis(:)),max(y_cross_axis(:))]);
        ylim([min(x_along_axis(:)),max(x_along_axis(:))]);
        colorbar;
    end
    
    amp = 20*log10(abs(slcStatic));
    pha = angle(slcStatic);
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1);
    c=surf(y_cross_axis,x_along_axis,amp);
    view(0,90);
    axis equal
    set(c,'EdgeColor','none');
    title('Amplitude (Static)');
    xlim([min(y_cross_axis(:)),max(y_cross_axis(:))]);
    ylim([min(x_along_axis(:)),max(x_along_axis(:))]);
    colorbar;
    subplot(2,1,2);
    c=surf(y_cross_axis,x_along_axis,pha);
    view(0,90);
    axis equal
    set(c,'EdgeColor','none');
    title('Phase (Static)');
    xlim([min(y_cross_axis(:)),max(y_cross_axis(:))]);
    ylim([min(x_along_axis(:)),max(x_along_axis(:))]);
    colorbar;
end

end

