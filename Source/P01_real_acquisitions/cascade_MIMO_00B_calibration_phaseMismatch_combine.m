%% cascade_MIMO_00B_calbration_antenna_combine.m
%
% The following function loads a set of calibration MAT-file acquired after 
% processing acquisitions of a single TIDEP-01012 sensor. The relevant
% matrices are averaged and stored as a new MAT-file at the same location.
%
% -------------------------------------------------------------------------
%
% Input:
%   
% - path2calibMatDir: A string specifying the directory where a set of
%                   calibration MAT-files is stored.
%
%
% -------------------------------------------------------------------------
% by Andreas Baumann-Ouyang, ETH Zürich (03nd August 2022)

function [path2calibMatComb] = cascade_MIMO_00B_calibration_phaseMismatch_combine(path2calibMatDir)

plt_fig = 1; % Plotting Overview - Slow if activated.

calibMat = dir(fullfile(path2calibMatDir,'*.mat')); % Get all Calibration Matrices

if plt_fig
    f1 = figure('units','normalized','outerposition',[0 0 1 1],'Name','RX Mismatch');
    f2 = figure('units','normalized','outerposition',[0 0 1 1],'Name','TX Mismatch');
    colormap(parula);
    
end

RxMismatch_Set = [];
TxMismatch_Set = [];
AngleMat_Set =  [];
RangeMat_Set =  [];
PeakValMat_Set =  [];
Rx_fft_Set = [];

N_calibMat = length(calibMat);
col = lines(N_calibMat);

for cal_i = 1:N_calibMat
    data = load(fullfile(calibMat(cal_i).folder,calibMat(cal_i).name));
    if cal_i == 1
        params = data.params;
        [m,n] = size(data.calibResult.TxMismatch);
    end

    RxMismatch = data.calibResult.RxMismatch;
    
    TxMismatch = data.calibResult.TxMismatch(:);
    
    RxMismatch_Set(:,end+1) = RxMismatch;
    TxMismatch_Set(:,end+1) = TxMismatch;
    AngleMat_Set(:,:,end+1) = data.calibResult.AngleMat;
    RangeMat_Set(:,:,end+1) = data.calibResult.RangeMat;
    PeakValMat_Set(:,:,end+1) = data.calibResult.PeakValMat;
    Rx_fft_Set(:,:,:,end+1) = data.calibResult.Rx_fft;

    if plt_fig
        figure(f1); hold on
        plot(RxMismatch,'Color',col(cal_i,:));
        
        figure(f2); hold on
        plot(TxMismatch,'Color',col(cal_i,:));
    end

end

if plt_fig
    figure(f1);
    title(sprintf('TIDEP-01012 Sensor - RX Mismatch'));
    xlabel('RXA No.');
    ylabel('Mismatch [degree]');

    figure(f2);
    title(sprintf('TIDEP-01012 Sensor - TX Mismatch'));
    xlabel('TXA No.');
    ylabel('Mismatch [degree]');
end

% New Calibration Values with Median Values for the series
calibResult.AngleMat = median(AngleMat_Set,3);
calibResult.RangeMat = median(RangeMat_Set,3);
calibResult.PeakValMat = median(PeakValMat_Set,3);
calibResult.RxMismatch = median(RxMismatch_Set,2);
calibResult.TxMismatch = reshape(median(TxMismatch_Set,2),m,n);
calibResult.Rx_fft = median(Rx_fft_Set,4);

if plt_fig
    figure(f1);
    plot(calibResult.RxMismatch,'Color','k','LineWidth',3);
    box on 
    xlim([1,max(xlim)])
    figure(f2);
    plot(calibResult.TxMismatch(:),'Color','k','LineWidth',3);
    box on 
    xlim([1,max(xlim)])
end

path2calibMatComb = fullfile(path2calibMatDir,'phaseMismatchCalibration.mat');
save(path2calibMatComb, 'calibResult','params');

exportgraphics(f1,fullfile(path2calibMatDir,'RxMismatch.png'),'Resolution',600);
savefig(f1,fullfile(path2calibMatDir,'RxMismatch.fig'));

exportgraphics(f2,fullfile(path2calibMatDir,'TxMismatch.png'),'Resolution',600);
savefig(f2,fullfile(path2calibMatDir,'TxMismatch.fig'));

end


