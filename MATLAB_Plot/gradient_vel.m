nx = 101;
nz = 101;
dx = 0.01;
dz = 0.01;
x = [0:dx:(nx-1)*dx];
z = [0:dz:(nz-1)*dz];

fid=fopen(['./gradient/gradient_vel.bin'],'r');
vel=fread(fid,[nz nx],'float');
fclose(fid);

figure('units','normalized','position',[0,0,0.4,0.5])
imagesc(x, z, vel);
colormap('jet');

c = colorbar;
c.Label.String = 'Velocity, m/s';
c.Label.FontName='Times New Roman';
c.Label.FontSize=30;
c.FontSize = 30; 
c.FontName = 'Times New Roman';
% set(gca, 'XAxisLocation', 'top');
set(gca, 'FontSize', 30, 'Fontname', 'Times New Roman');
set(gca,'Units','normalized','Position',[0.1 0.17 0.7 0.8]);
% ax=gca;
% ax.XTick = [160,330,480,640,800]; 
% ax.YTick = [80,160];
% xtic=ax.XTick*6.30;
% ytic=ax.YTick*6.30;
% ax.XTickLabel=xtic;
% ax.YTickLabel=ytic;
% text(0.5, 0, '(a)', 'FontSize', 30, 'FontName', 'Times New Roman', ...
%      'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
%      'Units', 'normalized', 'Parent', gca);

xlabel('Distance, km', 'FontSize', 30, 'FontName', 'Times New Roman');
ylabel('Depth, km', 'FontSize', 30, 'FontName', 'Times New Roman');

% print("./picture/gradientvel", '-dpng', '-r500');
