function [aoi, desc] = set_aoi_scXobs(create_aoi,...
                                               cplx,...
                                               coher,...
                                               coher_evol,...
                                               y_axis,...
                                               x_axis,...
                                               aoi_is_circle)
if ~exist('aoi_is_circle','var')
    aoi_is_circle = zeros(create_aoi);
end

N_limit = 10000;
if size(cplx,2)>N_limit
    N_skip = ceil(size(cplx,2)/N_limit);
    cplx = cplx(:,1:N_skip);
end

%complex_s0 = cplx(:,1);
ampl_ori = (mean(abs(cplx),2));
filt_perc = prctile(ampl_ori(:),85);
idx_perc = find(ampl_ori==min(ampl_ori(ampl_ori>=filt_perc)));
ampl_log10 = log10(ampl_ori); % Logarithmic scale of dataset
cb_limits = [ampl_log10(idx_perc),max(ampl_log10(:))];
filter_ampl = ampl_log10<min(cb_limits);
ampl_log10(filter_ampl)=min(cb_limits);

aoi = cell(create_aoi,1);
desc = cell(create_aoi,1); 

for i= 1:create_aoi
    figure( 'name','Area of Interests',...
            'units','normalized',...
            'outerposition',[0 0 1 1]);
    
    % Axes 1
    ax(1) = subplot(2,2,1);
    xlimits = [min(y_axis(:)),max(y_axis(:))];
    ylimits = [min(x_axis(:)),max(x_axis(:))];
    data_vis = log10(abs(max(cplx,[],2)));%;complex_s0));
    text_vis = "Maximum Amplitude (log10(A))";
    plot_polar_range_azimuth_2D_AB_preAX(y_axis, x_axis, data_vis,xlimits,ylimits,'scatter');
    title(text_vis)
    hold on
    colorbar;
    climits_val = [prctile(data_vis,50),prctile(data_vis,98)];
    caxis(climits_val); set(gca,'CLim',climits_val);
    ax(1).Color=[0,0,0];

    % Axes 3
    ax(2) = subplot(2,2,3);
    plot_polar_range_azimuth_2D_AB_preAX(y_axis, x_axis, coher, xlimits,ylimits,'scatter');
    colorbar; climits_val = [prctile(coher,70),1]; 
    caxis(climits_val); set(gca,'CLim',climits_val);
    title("Coherence between acquisitions");
    ax(2).Color=[0,0,0];

    % Axes 2
    ax(3) = subplot(2,2,2);   
    plot_polar_range_azimuth_2D_AB_preAX(y_axis, x_axis, ampl_log10,xlimits,ylimits,'scatter');
    title("Amplitude image. (Filtered A>85%)")
    %title("Amplitude image. (A for N=Avg)")
    ax(3).Color=[0,0,0];
    colorbar;

    % Axes 4
    ax(4) = subplot(2,2,4);
    plot_polar_range_azimuth_2D_AB_preAX(y_axis, x_axis, coher_evol,xlimits,ylimits,'scatter');
    title(sprintf("Change of Coherence based on filtred data (coh >= %.1f)",min(climits_val)));
    caxis([-2,1]);
    colorbar('Ticks',[-2,-1,0,1],...
     'TickLabels',{'Higher Coherence at Start','High Coherence','Low Coherence','Higher Coherence at End'});
    set(ax(4),'CLim',[-2,1]);
    ax(4).Color=[0,0,0];
    linkaxes([ax(1),ax(2),ax(3),ax(4)],'xy')
    if exist("xlimits","var") && exist("ylimits","var")
        xlim(xlimits); ylim(ylimits);
    end
    if ax(1).YLim(1)>=1000
        for axi=1:4
            ax(axi).YAxis.Exponent = 0;
            ax(axi).XAxis.Exponent = 0;
            ax(axi).YAxis.TickLabelFormat = '%.0f';
            ax(axi).XAxis.TickLabelFormat = '%.0f';
            ax(axi).XTickLabelRotation = 45;
        end
    end
    

    subplot(2,2,2);
    if aoi_is_circle(i)
        [aoi{i}, desc{i}] = pt2circle();
    else
        aoi_i = drawpolygon();
        aoi{i} = aoi_i.Position;
        answer = inputdlg('Description of the selected bins:');
        desc{i} = answer{1};
    end
end

aoi = round_aoi(aoi,3);

end

