nx = 501;
nz = 140;
dx = 10;
dz = 10;
x = [0:dx:(nx-1)*dx];
z = [0:dz:(nz-1)*dz];

fid=fopen(['/media/lqq/Data/wang/codes/Eikonal_FSM/wang/mar_140_501_h15s'],'r');
vel=fread(fid,[nz nx],'float');
fclose(fid);

fid = fopen('./mar/traveltime_ref_mar.dat', 'r');
time0 = fread(fid, [nz nx], 'float');
fclose(fid);
fid = fopen('./mar/traveltime_FSM.dat', 'r');
time1 = fread(fid, [nz nx], 'float');
fclose(fid);
fid = fopen('./mar/traveltime_stagered_FSM.dat', 'r');
time2 = fread(fid, [nz nx], 'float');
fclose(fid);

level = [0.1: 0.1: 0.9];
figure('units','normalized','position',[0,0,0.6,0.5])
imagesc(x, z, vel);
hold on;

contour(x,z,time0,level,'k-','ShowText','off','LineWidth',1);
contour(x,z,time1,level,'g--','ShowText','off','LineWidth',1);
contour(x,z,time2,level,'r--','ShowText','off','LineWidth',1);
colormap('jet');

c = colorbar;
c.Label.String = 'velocity(m/s)';
c.Label.FontName='Times New Roman';
c.Label.FontSize=20;
c.FontSize = 20; 
c.FontName = 'Times New Roman';
% set(gca, 'XAxisLocation', 'top');
set(gca, 'FontSize', 20, 'Fontname', 'Times New Roman');
% set(gca,'Units','normalized','Position',[0.07 0.05 0.83 0.79]);
% ax=gca;
% ax.XTick = [160,320,480,640,800]; 
% ax.YTick = [80,160];
% xtic=ax.XTick*6.25;
% ytic=ax.YTick*6.25;
% ax.XTickLabel=xtic;
% ax.YTickLabel=ytic;
% text(0.5, 0, '(a)', 'FontSize', 20, 'FontName', 'Times New Roman', ...
%      'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
%      'Units', 'normalized', 'Parent', gca);

xlabel('Distance(m)', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Depth(m)', 'FontSize', 20, 'FontName', 'Times New Roman');

print("./picture/marmousi_error", '-dpng', '-r500');
