function [] = video_tracking_v2(path2video,name_proc,isIndirect, azimuth, time_T0)

run_debugging = 0;
fontSize = 30;
roi_type = 'polygon';

[path2main,~,~] = fileparts(path2video);

[~,name_aoi_circle,~] = fileparts(path2video);
path2main_out = fullfile(path2main,name_aoi_circle);
if ~isfolder(path2main_out)
    mkdir(path2main_out)
end
path_aoi_circle = fullfile(path2main_out,sprintf('%s_aoi_circle.mat',name_proc));
path_aoi_track = fullfile(path2main_out,sprintf('%s_aoi_track.mat',name_proc));
path_aoi_scale = fullfile(path2main_out,sprintf('%s_aoi_scale.mat',name_proc));
path_displ = fullfile(path2main_out,sprintf('%s_displ_%s.mat',name_proc, datestr(datetime,'yyyymmdd_HHMMSS')));

videoReader = VideoReader(path2video);
frameRate = videoReader.FrameRate;
frameNum = videoReader.NumFrames;

sSize = get(0,'screensize');
    
%% ------------------------------------------------------------------------
% Get Area of Interest for Feature Tracking
step_i = 1;
fprintf('Step % 2d: Get Area of Interest for Feature Tracking\n',step_i);
frame = readFrame(videoReader);
frame_rot = imrotate(frame,azimuth,'bilinear');
[X_local_dim, Y_local_dim,~] = size(frame_rot);
[XimgRot,YimgRot] = meshgrid(X_local_dim-1:-1:0,Y_local_dim-1:-1:0); % Inverse
YLimits = [min(YimgRot(:)),max(YimgRot(:))];
XLimits = [min(XimgRot(:)),max(XimgRot(:))];
fig = figure('units','normalized','outerposition',[0 0 1 1]);
image(YimgRot(:), XimgRot(:), frame_rot);
set(gca,'YDir','normal','YDir','normal');
axis equal

if strcmp(roi_type,'rectangle')
    bx = msgbox('Set lower left to upper right points of a rectangle (= 2 Points)');
    pause(1);
    object_aoi = drawrectangle;
    object_name = 'Rectangle';
    objectShape = round(object_aoi.Position);
    objectShape2 = objectShape;
elseif strcmp(roi_type,'circle') || strcmp(roi_type,'polygon')
    if ~isfile(path_aoi_track)
        if strcmp(roi_type,'circle')
            bx = msgbox(sprintf('Set points on circle (>3 Points): %s',name_proc));
        else
            bx = msgbox(sprintf('Set points of polygon (>3 Points): %s',name_proc));
        end
        waitfor(bx)
        object_aoi = drawpolyline;
        C = object_aoi.Position;
        if strcmp(roi_type,'circle')
            a=[C(:,1) C(:,2) ones(size(C(:,1)))] \ ...
              [-(C(:,1).^2+C(:,2).^2)];
            xc = -.5*a(1);
            yc = -.5*a(2);
            Radi  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
            objectShape = round([xc, yc, R]);
            objectShape2 = round([xc-Radi, yc-Radi, 2*Radi, 2*Radi]);
            object_name = 'Circle';
            save(path_aoi_track,...
                    'objectShape',...
                    'objectShape2',...
                    'object_name',...
                    'xc',...
                    'yc',...
                    'Radi');
        else
            objectShape = round(C);
            objectShape2 = [min(objectShape(:,1)),min(objectShape(:,2)),...
                            max(objectShape(:,1))-min(objectShape(:,1)),...
                            max(objectShape(:,2))-min(objectShape(:,2))];
            object_name = 'Polygon';
            save(path_aoi_track,...
                    'objectShape',...
                    'objectShape2',...
                    'object_name');
        end
    else
        load(path_aoi_track);
    end
else
    error('Unknown type of shape.');
end
if exist('bx','var')
    if isvalid(bx)
        close(bx);
    else
        clearvars bx
    end
end

patch(objectShape(:,1),objectShape(:,2),...
      'red',...
      'EdgeColor','y', ...
      'FaceColor','none',...
      'LineStyle','--',...
      'LineWidth',3);

close(fig);

%% ------------------------------------------------------------------------
% Get trackable features and plot them
step_i = step_i + 1;
fprintf('Step % 2d: Visualised Trackable Features\n',step_i);
fig = figure('units','normalized','outerposition',[0 0 1 1]);
hold on

image(YimgRot(:), XimgRot(:), frame_rot);

gray = im2gray(rot90(frame_rot,2));
points = detectMinEigenFeatures(gray,'ROI',objectShape2);

if strcmp(roi_type,'circle')
    idx = ((points.Location(:,1)-xc).^2+...
           (points.Location(:,2)-yc).^2<=Radi^2);
    points = cornerPoints(points.Location(idx,:));
elseif strcmp(roi_type,'polygon')
    idx = inpolygon(points.Location(:,1),points.Location(:,2),...
                    objectShape(:,1),objectShape(:,2));
    points = cornerPoints(points.Location(idx,:));
end

point_tracked.coordinate_y = points.Location(:,1);
point_tracked.coordinate_x = points.Location(:,2);
point_tracked.coordinate_valid = ones(size(points.Location(:,1)));

scatter(point_tracked.coordinate_y,point_tracked.coordinate_x,...
        '+','Color','y');

if strcmp(roi_type,'circle')
    drawcircle('Center',[xc,yc],...
               'Radius',Radi,...
               'MarkerSize',0.01,...
               'Color',[1,0.2,0.2],...
               'LineWidth',1,...
               'FaceAlpha',0);
elseif strcmp(roi_type,'polygon')
    drawpolygon('Position',objectShape,...
               'Color',[1,0.2,0.2],...
               'LineWidth',1,...
               'FaceAlpha',0);
end

title('Detected interest points');

%% ------------------------------------------------------------------------
% Get scale of image objects
step_i = step_i + 1;
fprintf('Step % 2d: Get Scale Factor\n',step_i);
if ~isfile(path_aoi_scale)
    mb = msgbox('1.) Draw a line with known metric length.');
    waitfor(mb);
    scale_line.pixel_coord = drawline;
    scale_line.pixel = sqrt(sum(diff(scale_line.pixel_coord.Position,[],1).^2));
    scale_line.meter = inputdlg('2.) How long is the defined line in [meter]?');
    if ~isnan(str2double(scale_line.meter{1}))
        scale_line.meter = str2double(scale_line.meter{1});
    end
    scale_line.pixel2meter = scale_line.meter / scale_line.pixel;

    save(path_aoi_scale,'scale_line');
else
    load(path_aoi_scale);
end

%% ------------------------------------------------------------------------
% Scale of real world
m2mm = 1000; % [m] to [mm]
pixel2mm =  scale_line.pixel2meter * m2mm;

YimgRotmm = YimgRot * pixel2mm;
XimgRotmm = XimgRot * pixel2mm;

step_i = step_i + 1;
fprintf('Step % 2d: Scale [pixel] to [mm] with %.3f [mm/px]\n',...
        step_i,pixel2mm);

%% ------------------------------------------------------------------------
% Get Circle if applicable
if isIndirect
    fig2 = figure('units','normalized','outerposition',[0 0 1 1]);

    image(YimgRotmm(:),XimgRotmm(:),frame_rot);
    set(gca,'YDir','normal', 'XDir','normal','FontSize',fontSize);
    xlabel('Y^{G} [mm]');
    ylabel('X^{G} [mm]');
    hold on
    axis equal

    if ~isfile(path_aoi_circle)
        bx = msgbox('Set points of polygon (>3 Points) on circle for indirect points');
        waitfor(bx);
        object_aoi = drawpolyline;
        C = object_aoi.Position;
        a=[C(:,1) C(:,2) ones(size(C(:,1)))] \ ...
                  [-(C(:,1).^2+C(:,2).^2)];
        yc = -.5*a(1);
        xc = -.5*a(2);
        Radi  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
       
        save(path_aoi_circle,'C',"xc",'yc','Radi');

        drawcircle('Center',[yc,xc],...
                   'Radius',Radi,...
                   'MarkerSize',0.01,...
                   'Color',[1,0.2,0.2],...
                   'LineWidth',1,...
                   'FaceAlpha',0);
        scatter(yc,xc,'xy');
    end

    close(fig2)

    set(0,'CurrentFigure',fig);
end

%% ------------------------------------------------------------------------
% Tracking points over all frames
step_i = step_i + 1;
fprintf('Step % 2d: Tracking Points\n',step_i);
tracker = vision.PointTracker('MaxBidirectionalError',1);
initialize(tracker,points.Location,rot90(frame_rot,2));
dN = 1;
for f_i = 1:dN:frameNum
    frame = read(videoReader,f_i);
    frame_rot = imrotate(frame,azimuth,'bilinear');
    frame_rot2 = rot90(frame_rot,2);
    [points,validity] = tracker(frame_rot2);
    
    point_tracked.coordinate_y(:,end+1) = points(:,1);
    point_tracked.coordinate_x(:,end+1) = points(:,2);
    point_tracked.coordinate_valid(:,end+1) = validity;
    if run_debugging
        fig3 = figure;
        image(YimgRot(:), XimgRot(:), frame_rot); hold on
        scatter(point_tracked.coordinate_y(:,end),...
                point_tracked.coordinate_x(:,end),...
                '.y');
        set(0,'CurrentFigure',fig)
    end

end

idx = (frameNum+1)-sum(point_tracked.coordinate_valid,2);
idx = idx==0;

point_tracked.coordinate_y =  point_tracked.coordinate_y(idx,:);
point_tracked.coordinate_x =  point_tracked.coordinate_x(idx,:);
point_tracked.coordinate_valid =  point_tracked.coordinate_valid(idx,:);
point_tracked.coordinate_time = time_T0 + seconds(1/frameRate * [0:frameNum]);

scatter(point_tracked.coordinate_y(:,1),point_tracked.coordinate_x(:,1),...
        'o','Color','y');

set(gca,'FontSize',fontSize);
xlabel('Y^G [pixel]');
ylabel('X^G [pixel]');
axis equal
xlim(YLimits);
ylim(XLimits);

exportgraphics(fig,replace(path_aoi_track,'.mat','.png'),'Resolution',600);
saveas(fig,replace(path_aoi_track,'.mat','.fig'));

close(fig);

%% ------------------------------------------------------------------------
% Scale to Real World
point_tracked.coordinate_y = (point_tracked.coordinate_y) * pixel2mm;
point_tracked.coordinate_x = (point_tracked.coordinate_x) * pixel2mm;


%% ------------------------------------------------------------------------
% Figure of 2D-Displacement Trajectories
step_i = step_i + 1;
fprintf('Step % 2d: Visualise Displacement Trajectories\n',step_i);
fig_traj = figure('units','normalized','outerposition',[0 0 1 1]);
videoReader = VideoReader(path2video);
frame = readFrame(videoReader);
frame_rot = imrotate(frame,azimuth,'bilinear');

image(YimgRotmm(:),XimgRotmm(:),frame_rot);
set(gca,'YDir','normal', 'XDir','normal','FontSize',fontSize);
xlabel('Y^{G} [mm]');
ylabel('X^{G} [mm]');

hold on
scatter(point_tracked.coordinate_y(:,1),point_tracked.coordinate_x(:,1),...
        80,...
        '+y',...
        'DisplayName','Tracked Features');
avg_point_y = mean(point_tracked.coordinate_y,1);
avg_point_x = mean(point_tracked.coordinate_x,1);
plot(avg_point_y,avg_point_x,'g','LineWidth',2,...
        'DisplayName','Trajectory of Averaged Feature Locations');
scatter(avg_point_y(1),avg_point_x(1),100,'^r','filled',...
        'DisplayName','Averaged Feature Locations at Start');
scatter(avg_point_y(end),avg_point_x(end),100,'vr','filled',...
        'DisplayName','Averaged Feature Locations at End');
grid on
box on 
axis equal
%l = legend('AutoUpdate','off','Location','northoutside','NumColumns',3);
y_lim = [min(YimgRotmm(:)),max(YimgRotmm(:))];
x_lim = [min(XimgRotmm(:)),max(XimgRotmm(:))];
xlim(y_lim);
ylim(x_lim);

% Get Coordinates for North-Symbol
north_arrow(1,1) = y_lim(1)+diff(y_lim)*0.05;
north_arrow(1,2) = x_lim(1)+diff(x_lim)*0.05;
north_arrow(2,:) = [0; diff(x_lim)*0.05];
arrow_length = min( [sqrt(sum(north_arrow(2,:).^2)),...
                    sqrt( (north_arrow(1,1)-north_arrow(1,1)+north_arrow(2,1)).^2 + ...
                          (north_arrow(1,2)-north_arrow(1,2)+north_arrow(2,2)).^2)*0.99]);

try
    davinci('arrow',...
            'X', [north_arrow(1,1), north_arrow(1,1)+north_arrow(2,1)],...
            'Y', [north_arrow(1,2), north_arrow(1,2)+north_arrow(2,2)],...
            'Head.Length', arrow_length,...
            'Head.Sweep', 0,...
            'Head.Width', arrow_length*0.75,...,...
            'Color','w');
catch
    fprintf('Please install davinci_draw_R2017a to visualise North-Arrow\n');
end
text([north_arrow(1,1)+north_arrow(2,1)*1.2],...
     [north_arrow(1,2)+north_arrow(2,2)*1.5],...
     'N',...
     'HorizontalAlignment','center',...
     'FontSize',fontSize,...
     'Color','w');

% Get Coordinates for inverse image axis
y0 = min(YimgRotmm(:));
x0 = min(XimgRotmm(:));
y1 = max(YimgRotmm(:));
x1 = max(XimgRotmm(:));
frame_rot_tmp = sum(frame_rot,3);
frame_rot_tmp(frame_rot_tmp==0)=nan;
x2 = frame_rot_tmp(:,1);
x2 = mean(find(x2>0)) * pixel2mm; % X^G bottom left
x3 = frame_rot_tmp(:,end);
x3 = mean(find(x3>0)) * pixel2mm; % X^G top right
y2 = frame_rot_tmp(1,:);
y2 = mean(find(y2>0)) * pixel2mm; % Y^G top left
y3 = frame_rot_tmp(end,:);
y3 = mean(find(y3>0)) * pixel2mm; % Y^G top left

C_bl = [y0,x2]; % Coordinate Bottom-Left
C_tr = [y1,x3]; % Coordinate Top-Right
C_br = [y2,x0]; % Coordinate Bottom-Right
C_tl = [y3,x1]; % Coordinate Top-Left
dax = (y1-y0)*0.1; % length of axis

north_arrow(1,:) = C_bl;
Rot_matArrow = [cosd(azimuth) sind(azimuth);...
          -sind(azimuth) cosd(azimuth)]; % Clockwise
north_arrow(2,:) = Rot_matArrow * [0; dax];
try
    davinci('arrow',...
            'X', [north_arrow(1,1), north_arrow(1,1)+north_arrow(2,1)*1.02],...
            'Y', [north_arrow(1,2), north_arrow(1,2)+north_arrow(2,2)*1.02],...
            'Head.Length', dax*0.1*1.5,...
            'Head.Width',dax*0.1*1.8,...
            'Head.Sweep', 0,...
            'Shaft.Width',dax*0.035,...
            'Color','k');
    davinci('arrow',...
            'X', [north_arrow(1,1), north_arrow(1,1)+north_arrow(2,1)],...
            'Y', [north_arrow(1,2), north_arrow(1,2)+north_arrow(2,2)],...
            'Head.Length', dax*0.1,...
            'Head.Width',dax*0.1,...
            'Head.Sweep', 0,...
            'Shaft.Width',dax*0.025,...
            'Color','w');
catch
    fprintf('Please install davinci_draw_R2017a to visualise North-Arrow\n');
end

text([north_arrow(1,1)+north_arrow(2,1)*1],...
     [north_arrow(1,2)+north_arrow(2,2)*1.35],...
     'X^{Ci}',...
     'HorizontalAlignment','center',...
     'FontSize',fontSize,...
     'FontWeight','bold',...
     'Color','k');
text([north_arrow(1,1)+north_arrow(2,1)*1],...
     [north_arrow(1,2)+north_arrow(2,2)*1.35],...
     'X^{Ci}',...
     'HorizontalAlignment','center',...
     'FontSize',fontSize,...
     'Color','w');

Rot_matArrow = [cosd(90+azimuth) sind(90+azimuth);...
          -sind(90+azimuth) cosd(90+azimuth)]; % Clockwise
north_arrow(2,:) = Rot_matArrow * [0; dax];
try
    davinci('arrow',...
            'X', [north_arrow(1,1), north_arrow(1,1)+north_arrow(2,1)*1.02],...
            'Y', [north_arrow(1,2), north_arrow(1,2)+north_arrow(2,2)*1.02],...
            'Head.Length', dax*0.1*1.5,...
            'Head.Width',dax*0.1*1.8,...
            'Head.Sweep', 0,...
            'Shaft.Width',dax*0.035,...
            'Color','k');
    davinci('arrow',...
            'X', [north_arrow(1,1), north_arrow(1,1)+north_arrow(2,1)],...
            'Y', [north_arrow(1,2), north_arrow(1,2)+north_arrow(2,2)],...
            'Head.Length', dax*0.1,...
            'Head.Width',dax*0.1,...
            'Head.Sweep', 0,...
            'Shaft.Width',dax*0.025,...
            'Color','w');
catch
    fprintf('Please install davinci_draw_R2017a to visualise North-Arrow\n');
end
text([north_arrow(1,1)+north_arrow(2,1)*1.25],...
     [north_arrow(1,2)+north_arrow(2,2)*2.25],...
     'Y^{Ci}',...
     'HorizontalAlignment','center',...
     'FontSize',fontSize,...
     'FontWeight','bold',...
     'Color','k');
text([north_arrow(1,1)+north_arrow(2,1)*1.25],...
     [north_arrow(1,2)+north_arrow(2,2)*2.25],...
     'Y^{Ci}',...
     'HorizontalAlignment','center',...
     'FontSize',fontSize,...
     'Color','w');

title('Tracked Feature Points');

exportgraphics(fig_traj,replace(path_displ,'.mat','_Traj.png'),'Resolution',600);
saveas(fig_traj,replace(path_displ,'.mat','_Traj.fig'));

%% ------------------------------------------------------------------------
% Figure of East-West and North-South Components
step_i = step_i + 1;
fprintf('Step % 2d: Display X^G and Y^G components \n',step_i);
fig = figure('units','normalized','outerposition',[0 0 1 1]);

avg_point_y_rot = avg_point_y - mean(avg_point_y(1:10)); % Up-Down
avg_point_x_rot = avg_point_x - mean(avg_point_x(1:10)); % Left-Right

axmn(1) = subplot(2,1,1); hold on
T = seconds(point_tracked.coordinate_time - min(point_tracked.coordinate_time));
yline(0,'LineWidth',1.5,'Color','k');
plot(T,avg_point_y_rot,'LineWidth',2);
ylabel('\DeltaD_{X^G} [mm]');
ylim_val1 = ylim();
box on

axmn(2) = subplot(2,1,2); hold on
yline(0,'LineWidth',1.5,'Color','k');
plot(T,avg_point_x_rot,'LineWidth',2);
ylabel('\DeltaD_{Y^G} [mm]');
xlabel('Time [s]');
ylim_val2 = ylim();
box on

ylim_valmax = max(abs([ylim_val1,ylim_val2]));
for axi=1:2
    set(axmn(axi),'YLim',[-ylim_valmax,ylim_valmax],...
                  'XLim',[min(T),...
                          max(T)],...
                  'FontSize',fontSize);

end
sgtitle('Movement of Averaged Features','FontSize',fontSize*1.1);

exportgraphics(fig(end),replace(path_displ,'.mat','_AvgCoordComp.png'),...
               'Resolution',600);
saveas(fig(end),replace(path_displ,'.mat','_AvgCoordComp.fig'));

close(fig)

%% ------------------------------------------------------------------------
% Points on Circular Object/Tower
if isIndirect
    frame = read(videoReader,1);
    frame_rot = imrotate(frame,azimuth,'bilinear');
    fig = figure('units','normalized','outerposition',[0 0 1 1]);

    image(YimgRotmm(:),XimgRotmm(:),frame_rot);
    set(gca,'YDir','normal', 'XDir','normal','FontSize',fontSize);
    xlabel('Y^{G} [mm]');
    ylabel('X^{G} [mm]');
    hold on
    axis equal
    xlim(y_lim);
    ylim(x_lim);

    if isfile(path_aoi_circle)
        load(path_aoi_circle)
    else
        bx = msgbox('Set points of polygon (>3 Points)');
        waitfor(bx);
        object_aoi = drawpolyline;
        C = object_aoi.Position;
        a=[C(:,1) C(:,2) ones(size(C(:,1)))] \ ...
                  [-(C(:,1).^2+C(:,2).^2)];
        yc = -.5*a(1);
        xc = -.5*a(2);
        Radi  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
        save(path_aoi_circle,'C',"xc",'yc','Radi');
    end
    drawcircle('Center',[yc,xc],'Radius',Radi,'MarkerSize',0.01,'Color',[1,0.2,0.2],'LineWidth',1,'FaceAlpha',0);
    
    scatter(yc,xc,'xy');

    

    % Calculate Transformations
    n_Pt = size(point_tracked.coordinate_y,1);
    n_Time = size(point_tracked.coordinate_y,2);
    n_coord = 2;
    
    Tc = zeros(n_Time,size(C,1)+1,n_coord+1);
    Tc(1,1:end-1,1:2) = C;
    Tc(1,end,1:2) = [yc,xc];
    Tc(1,:,3) = ones(size(C,1)+1,1);
    
    tower_center.coordinate_y = zeros(n_Time,1);
    tower_center.coordinate_y(1) = yc;
    tower_center.coordinate_x = zeros(n_Time,1);
    tower_center.coordinate_x(1) = xc;

    tower_center.radi = zeros(n_Time,1);
    tower_center.radi(1) = Radi;
    
    for ti = 1:size(point_tracked.coordinate_y,2)
        idx = 1;
        T0 = double([point_tracked.coordinate_y(:,idx),...
                     point_tracked.coordinate_x(:,idx)]);
        T1 = double([point_tracked.coordinate_y(:,ti),...
                     point_tracked.coordinate_x(:,ti)]);
    
        T0 = [T0,ones(n_Pt,1)]';
        T1 = [T1,ones(n_Pt,1)]';
    
        A = round(T1/T0,8);
    
        Tc_pre = squeeze(Tc(idx,:,:));
        Tc(ti,:,:) = [A*Tc_pre']';
        
        C = squeeze(Tc(ti,1:end-1,1:2)); % Last one is center coordinate!
        a=[C(:,1) C(:,2) ones(size(C(:,1)))] \ ...
              [-(C(:,1).^2+C(:,2).^2)];
    
        tower_center.coordinate_y(ti) = -.5*a(1);
        tower_center.coordinate_x(ti) = -.5*a(2);
        tower_center.radi(ti) =  sqrt((a(1)^2+a(2)^2)/4-a(3));
        if ti==2
            tower_center.coordinate_time = point_tracked.coordinate_time';
        end       
    end

    scatter(tower_center.coordinate_y(1),...
            tower_center.coordinate_x(1),'oy');
    pause(5);
    close(fig);

    tower_center.coordinate_N = tower_center.coordinate_x;
    tower_center.coordinate_N = tower_center.coordinate_N - mean(tower_center.coordinate_N);
    tower_center.coordinate_E = tower_center.coordinate_y;
    tower_center.coordinate_E = tower_center.coordinate_E - mean(tower_center.coordinate_E);

    N_avg = 100;
    E_C = mean(tower_center.coordinate_y(1:N_avg));
    N_C = mean(tower_center.coordinate_x(1:N_avg));
    dNCoord = tower_center.coordinate_N-mean(tower_center.coordinate_N(1:N_avg));

    north_arrow(1,1) = y_lim(2)-diff(y_lim)*0.05-E_C;
    north_arrow(1,2) = x_lim(2)-diff(x_lim)*0.12-N_C;
    north_arrow(2,:) = [0; diff(x_lim)*0.05];
    arrow_length = min( [sqrt(sum(north_arrow(2,:).^2)),...
                        sqrt( (north_arrow(2,1)).^2 + ...
                              (north_arrow(2,2)).^2)*0.99]);

    step_i = step_i + 1;
    fprintf('Step % 2d: Create Video-Tracking-Animation (Tower Center)\n',step_i);
    XLim_Val = [min([tower_center.coordinate_y;YimgRotmm(:)]),...
                max([tower_center.coordinate_y;YimgRotmm(:)])]-E_C;
    YLim_Val = [min([tower_center.coordinate_x;XimgRotmm(:)]),...
                max([tower_center.coordinate_x;XimgRotmm(:)])]-N_C;


    fig_traj2 = figure('units','pixels',...
                        'outerposition',[0 0 sSize(3) sSize(4)],...
                        'visible','off');
    videoReader = VideoReader(path2video);
    videoWriter = VideoWriter(replace(path_displ,'.mat','_Tracking_02.avi'));
    videoWriter.FrameRate = videoReader.FrameRate;
    open(videoWriter);
    dN = floor(videoReader.FrameRate/4);
    step_i = step_i + 1;
    fprintf('Step % 2d: Start Framing...\n',step_i);
    time0 = min(tower_center.coordinate_time);
    Time_plot = seconds(tower_center.coordinate_time-time0);
    for f_i = 1:dN:frameNum

        set(0,'CurrentFigure',fig_traj2)

        frame = read(videoReader,f_i);
        frame_rot = imrotate(frame,azimuth,'bilinear');
        subplot(15,1,1:11);
        imagesc(YimgRotmm(:)-E_C,XimgRotmm(:)-N_C,frame_rot);
        set(gca,'YDir','normal',...
                'XDir','normal',...
                'Color','k',...
                'FontSize',fontSize);
        ylabel('\DeltaD_{X^{G}} [mm]');
        xlabel('\DeltaD_{Y^{G}} [mm]');
        axis equal
        xlim(XLim_Val);
        ylim(YLim_Val);
             
        if run_debugging
           xline(0,'w');
           yline(0,'w');
        end
        
        hold on
        C = squeeze(Tc(f_i,1:end-1,1:2));

        drawcircle('Center',[tower_center.coordinate_y(f_i)-E_C,...
                             tower_center.coordinate_x(f_i)-N_C],...
                    'Radius',tower_center.radi(f_i),...
                    'MarkerSize',0.01,...
                    'Color',[1,0.2,0.2],...
                    'LineWidth',1,...
                    'FaceAlpha',0);

        scatter(C(:,1)-E_C,C(:,2)-N_C,...
                'xy',...
                'DisplayName','Tracked Tower Circle');
        scatter(point_tracked.coordinate_y(:,f_i)-E_C,...
                point_tracked.coordinate_x(:,f_i)-N_C,...
                '.y',...
                'DisplayName','Tracked Turbine Features');
        scatter(tower_center.coordinate_y(f_i)-E_C,...
                tower_center.coordinate_x(f_i)-N_C,...
                '+y',...
                'DisplayName','Tracked Tower Center');
        scatter(tower_center.coordinate_y(f_i)-E_C,...
                tower_center.coordinate_x(f_i)-N_C,...
                'oy',...
                'DisplayName','Tracked Tower Center');
        
        try
            davinci('arrow',...
                'X', [north_arrow(1,1), north_arrow(1,1)+north_arrow(2,1)],...
                'Y', [north_arrow(1,2), north_arrow(1,2)+north_arrow(2,2)],...
                'Head.Length', arrow_length,...
                'Head.Sweep', 0,...
                'Head.Width', arrow_length*0.75,...,...
                'Color','w');
        catch
            fprintf('Please install davinci_draw_R2017a to visualise North-Arrow\n');
        end
        text(north_arrow(1,1)+north_arrow(2,1)*1.2,...
             north_arrow(1,2)+north_arrow(2,2)*1.5,...
             'N',...
             'HorizontalAlignment','center',...
             'FontSize',fontSize,...
             'Color','w');
       
        title(sprintf("Video Camera (%s)",datestr(tower_center.coordinate_time(f_i),'dd.mm.yyyy HH:MM:SS')))
        hold off

        subplot(15,1,14:15);

        plot(Time_plot, dNCoord,'k'); hold on
        scatter(Time_plot(f_i), dNCoord(f_i),'or','filled');
        ylabel('X^G [mm]');
        xlabel('Time [s]');
        xlim([min(Time_plot),...
              max(Time_plot)]);
        dyDist = [min(dNCoord),max(dNCoord)];
        dyDist = dyDist + diff(dyDist)*0.1 * [-1 1];
        ylim(dyDist)
        set(gca,'FontSize',fontSize);
        hold off
        xtime_ticks = xticks();

        set(fig_traj2,...
            'units','pixels',...
            'outerposition',[0 0 sSize(3) sSize(4)])

        frameW = getframe(fig_traj2);
        writeVideo(videoWriter, frameW); 

        if min(abs(Time_plot(f_i)-xtime_ticks))*videoReader.FrameRate<(dN/2)
            str_txt = sprintf('Video_%04ds',...
                             round(Time_plot(f_i)));
            exportgraphics(fig_traj2,...
                           replace(path_displ,'.mat',[str_txt,'.png']),...
                           'Resolution',600);
            saveas(fig_traj2,replace(path_displ,'.mat',[str_txt,'.fig']));
        end

        if f_i == 1
            % Update Trajectory Figure
            set(0,'CurrentFigure',fig_traj)

            drawcircle(...
                    'Center',[tower_center.coordinate_y(f_i),...
                              tower_center.coordinate_x(f_i)],...
                    'Radius',tower_center.radi(f_i),...
                    'MarkerSize',0.01,...
                    'Color','yellow',...
                    'LineWidth',2,...
                    'FaceAlpha',0,...
                    'HandleVisibility','off');

            scatter(C(:,1),...
                    C(:,2),...
                    100,...
                    'xy',...
                    'DisplayName','Tracked Tower Circle');
            scatter(C(:,1),C(:,2),100,...
                'oy',...
                'HandleVisibility','off');

            scatter(...
                    tower_center.coordinate_y(f_i),...
                    tower_center.coordinate_x(f_i),...
                    120,...
                    '+y',...
                    'HandleVisibility','off');

            yLim = ylim()-N_C;
            xLim = xlim()-E_C;
            yLim0 = 0:-2000:ceil(yLim(1));
            yLim1 = 0:2000:floor(yLim(2));
            yLimN = [flip(yLim0(2:end)),yLim1];
            yLimNV = yLimN+N_C;
            yticks(yLimNV);
            yticklabels(yLimN);
            xLim0 = 0:-2000:ceil(xLim(1));
            xLim1 = 0:2000:floor(xLim(2));
            xLimN = [flip(xLim0(2:end)),xLim1];
            xLimNV = xLimN+E_C;
            xticks(xLimNV);
            xticklabels(xLimN);

            xlabel('\DeltaD_{Y^G} [mm]');
            ylabel('\DeltaD_{X^G} [mm]');
            exportgraphics(fig_traj,replace(path_displ,'.mat','_Traj2.png'),'Resolution',600);
            saveas(fig_traj,replace(path_displ,'.mat','_Traj2.fig'));

            set(0,'CurrentFigure',fig_traj2)
        end
    end
    close(videoWriter);
    step_i = step_i + 1;
    fprintf('Step % 2d: ...End of Framing\n',step_i);
    
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,1,1);
    plot(point_tracked.coordinate_time,tower_center.coordinate_N);
    ylabel('X_{North} [mm]');
    xlabel('Time');
    subplot(2,1,2);
    plot(point_tracked.coordinate_time,tower_center.coordinate_E);
    ylabel('Y_{East} [mm]');
    
    exportgraphics(fig,replace(path_displ,'.mat','_CenterCoordComp.png'),'Resolution',600);
    saveas(fig,replace(path_displ,'.mat','_CenterCoordComp.fig'));
else
    tower_center.coordinate_N = avg_point_y_rot;
    tower_center.coordinate_N = tower_center.coordinate_N - mean(tower_center.coordinate_N);
    tower_center.coordinate_E = avg_point_x_rot;
    tower_center.coordinate_E = tower_center.coordinate_E - mean(tower_center.coordinate_E);
    tower_center.coordinate_time = point_tracked.coordinate_time';
end

save(path_displ,'tower_center');

end