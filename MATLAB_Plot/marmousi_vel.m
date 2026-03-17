nx = 501;
nz = 140;
dx = 0.01;
dz = 0.01;
x = [0:dx:(nx-1)*dx];
z = [0:dz:(nz-1)*dz];

fid=fopen(['./marmousi/mar_140_501_h15s'],'r');
vel=fread(fid,[nz nx],'float');
fclose(fid);

figure('units','normalized','position',[0.05,0.05,0.4,0.3])
imagesc(x, z, vel);
colormap('jet');

c = colorbar;
c.Label.String = 'Velocity, m/s';
c.Label.FontName='Times New Roman';
c.Label.FontSize=25;
c.FontSize = 25; 
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

xlabel('Distance, km', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Depth, km', 'FontSize', 20, 'FontName', 'Times New Roman');
set(gca, 'YDir', 'reverse', 'FontSize', 25);

% print("./picture/marmousivel", '-dpng', '-r500');
