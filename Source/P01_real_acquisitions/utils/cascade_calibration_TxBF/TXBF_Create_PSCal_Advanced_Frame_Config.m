%% TXBF_Create_PSCal_Advanced_Frame_Config.m
% script to convert the TX phase shifter calibration matrix data
% to a form that can be used by the existing TXBF Advanced Frame
% Configuration example setup. 
%
% phaseOffsetValues = [numDevices, numTXPerDevice, numRXPerStystem, numPSOffsets]
% Ph = [numPSOffsets - 1, numRXPerSystem, bum]
%

function [path2psCal] = TXBF_Create_PSCal_Advanced_Frame_Config(phaseShiftCalFile)

    load(phaseShiftCalFile,'*');
    
    if ~exist('phaseOffsetValues','var')
        error('MAT-File did not have a variable phaseOffsetValues');
    end

    numTX = 12;
    numRX = 16;
    numPSOffsets = 64;
    
    
    %% map phaseOffsetValues to Ph matrix format
    
    Ph = zeros(numPSOffsets - 1, numRX, numTX);
    
    for idxTX = 1:numTX
        for idxRX = 1:numRX
           for idxPSOffsets = 1:numPSOffsets - 1
    
               Ph(idxPSOffsets, idxRX, idxTX) = phaseOffsetValues(idxTX, idxRX, idxPSOffsets + 1);
    
           end 
        end
    end   
    
    
    %% graph both cal matrices
    fig(1) = figure('Name','Phase Offset Values',...
            'Units','centimeters',...
            'Position',[0,0,45,15]);  % plot phase offset values
    ax(1) = subplot(1,2,1);
    ax(2) = subplot(1,2,2);
    
    
    for idxTX = 1:numTX
        for idxRX = 1%:numRX * numDevices                
            plot(ax(1), squeeze(phaseOffsetValues(idxTX, idxRX, :)));
            hold(ax(1), 'on');
        end
    end   
    
    
    for idxTX = 1:numTX 
        for idxRX = 1%:numRX * numDevices
                plot(ax(2), squeeze(Ph(:, idxRX, idxTX)));
                hold(ax(2), 'on');
        end
    end   
    
    title(ax(1), 'phaseOffsetValues(idxDevices, idxTX, idxRX, idxPhaseOffsets)');
    title(ax(2), 'Ph(idxPhaseOffsets, idxRX, idxTX)');
    
    % save calibration data to .mat file
    path2psCal = fullfile(fileparts(phaseShiftCalFile), 'phaseShifterCalibration.mat');
    save(path2psCal, 'Ph');

    savefig(fig(1),replace(path2psCal,'.mat','.fig'));
    exportgraphics(fig(1),...
            replace(path2psCal,'.mat','.png'),...
            'Resolution',600);

end