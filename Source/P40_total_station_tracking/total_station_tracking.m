function [total_station] = total_station_tracking(path2txt,year, time_shift, azimuth)

[path2main,name2txt,~] = fileparts(path2txt);

csv_data = readtable(path2txt);
csv_fields = fields(csv_data);

for fi = 1:length(csv_fields)-3
    datastring = csv_data.(csv_fields{fi});

    if strcmp(csv_fields{fi},'Punktnummer')
        data = datastring(1:end-1);
        num_data = size(data,1);
        if iscell(data)
            data_tmp = repmat('A',size(data,1),12);
            for ci = 1:num_data
                data_tmp(ci,:) = data{ci}(5:end);
            end
            data = data_tmp;
        else

            data = num2str(data);
            data = data(:,5:end);
        end
        

        if length(year) == 1
            year = repmat(year,num_data,1);
        end

        month = str2num(num2str(data(:,1:2)));
        day = str2num(num2str(data(:,3:4)));
        hour = str2num(num2str(data(:,5:6)));
        minu = str2num(num2str(data(:,7:8)));
        sec = str2num(num2str(data(:,9:10)));
        millisec = str2num(num2str(data(:,11:12)));

        total_station.Time = datetime(year,...
                                      month,...
                                      day,...
                                      hour,...
                                      minu,...
                                      sec,...
                                      millisec);
        total_station.Time = total_station.Time - time_shift;

        total_station.ID = [1:num_data]';

        if total_station.Time(1)>total_station.Time(end)
            total_station.ID = sort(total_station.ID,'descend');
        end
    elseif strcmp(csv_fields{fi},'V_Winkel')
        data = split(datastring(1:end-1),' ');
        total_station.zenith = zeros(size(total_station.ID));
        for vi = 1:num_data
            total_station.zenith(vi) = str2num(data{vi});
        end
    elseif strcmp(csv_fields{fi},'Hz_Winkel')
        data = split(datastring(1:end-1),' ');
        total_station.azimuth = zeros(size(total_station.ID));
        for vi = 1:num_data
            total_station.azimuth(vi) = str2num(data{vi});
        end
    elseif strcmp(csv_fields{fi},'Schr_gdistanz')
        total_station.s_dist = datastring(1:end-1);
    end
end

total_station.s_diff = total_station.s_dist-mean(total_station.s_dist);
total_station.s_diff_hz = sind(total_station.zenith) .* total_station.s_diff;
total_station.s_diff_N = cosd(azimuth) * total_station.s_diff_hz;
total_station.s_diff_E = sind(azimuth) * total_station.s_diff_hz;

figure('Units','normalized','Position',[0,0,1,1]);

ax(1) = subplot(3,1,1);
hold on
grid on
box on
plot(total_station.Time,...
     total_station.s_diff);
xlabel('Time');
ylabel('\Delta D_{LOS} [m]');

ax(2) = subplot(3,1,2);
hold on
grid on
box on
plot(total_station.Time,...
     total_station.s_diff_N);
xlabel('Time');
ylabel('\Delta D_{X^G} [m]');

ax(3) = subplot(3,1,3);
hold on
grid on
box on
plot(total_station.Time,...
     total_station.s_diff_E);
yline(0,'LineWidth',1,'LineStyle','--');
xlabel('Time');
ylabel('\Delta D_{X^G} [m]');

for ax_i = 1:length(ax)
    yline(ax(ax_i),0,'LineWidth',1,'LineStyle','--');
    set(ax(ax_i),...
        'XLim',[min(total_station.Time),max(total_station.Time)],...
        'YLim',[min([total_station.s_diff(:);...
                     total_station.s_diff_N(:);...
                     total_station.s_diff_E(:)]),...
                max([total_station.s_diff(:);...
                     total_station.s_diff_N(:);...
                     total_station.s_diff_E])],...
        'LineWidth',2);
end

Link = linkprop(ax,{'XLim', 'YLim'});
setappdata(gcf, 'StoreTheLink', Link);

path2out = fullfile(path2main,name2txt);

if ~isfolder(path2out)
    mkdir(path2out)
end

pathname = fullfile(path2out,'total_station');
save(pathname,'total_station');
saveas(gcf,pathname);
exportgraphics(gcf,[pathname,'.png'],'Resolution',600);

end

