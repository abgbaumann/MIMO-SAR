function [aoi_c, desc_c] = pt2circle()

msgbox("First Input is the Center of the Circle. All following are points on the circle. The furthest one defines the radius.");
aoi_draw = drawpolyline();
aoi = aoi_draw.Position;

coord_center = aoi(1,:);

coord_radi = aoi(2:end,:);

dR = sqrt(sum((coord_radi-coord_center).^2,2));

dR = max(dR);
dradi = 1; % degree

radi = 0:dradi:360-dradi;
radi = radi*pi/180;

aoi_c(:,1) = coord_center(1) + cos(radi) * dR;
aoi_c(:,2) = coord_center(2) + sin(radi) * dR;

answer = inputdlg('Description of the selected bins:');
desc_c = answer{1};

end

