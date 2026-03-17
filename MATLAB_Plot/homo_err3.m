nx = 101;
nz = 101;
dx = 0.01;
dz = 0.01;
x = [0:dx:(nx-1)*dx];
z = [0:dz:(nz-1)*dz];

fid=fopen(['./homo/traveltime_homo_analy.bin'],'r');
ana=fread(fid,[nz nx],'float');
fclose(fid);

% fid=fopen(['/media/lqq/Data/wang/codes/Eikonal_FSM/wang/traveltime_FSM.dat'],'r');
% fid=fopen(['/media/lqq/Data/wang/codes/Eikonal_FSM/wang/traveltime_2ndFSM.dat'],'r');
fid=fopen(['./homo/traveltime_stagered_FSM.dat'],'r');
fsm=fread(fid,[nz nx],'float');
fclose(fid);

error_abs = abs(fsm-ana);
error_re = (error_abs ./ ana)*100; 
fprintf("%f %f", min(error_abs(:)), max(error_abs(:)));
levels = 0.0:0.5:2;

figure('units','normalized','position',[0,0,0.9,0.45]);

subplot('Position', [0.05 0.15 0.25 0.8]);
contour(x,z,ana,'k-','ShowText','off','LevelStep',0.1,'LineWidth',1);
hold on;
contour(x,z,fsm,'r--','ShowText','off','LevelStep',0.1,'LineWidth',1);
hold off;
xlabel('Distance, km', 'FontSize', 25, 'FontName', 'Times New Roman');
ylabel('Depth, km', 'FontSize', 25, 'FontName', 'Times New Roman');
text(-0.1, 1.05, '(g)', 'FontSize', 30, 'FontName', 'Times New Roman', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'Units', 'normalized', 'Parent', gca, 'FontWeight', 'bold');
set(gca, 'YDir', 'reverse', 'FontSize', 25);
legend('Analytical','RSGFSM', 'FontSize', 25, 'FontName', 'Times New Roman', Location='south',Orientation='horizontal');

subplot('Position', [0.36 0.15 0.25 0.8]);
imagesc(x, z, error_abs);
colormap('gray');
c = colorbar;
c.Label.String = 'Absolute error, s';
c.Label.FontName='Times New Roman';
c.Label.FontSize=25;
c.FontSize = 25; 
c.FontName = 'Times New Roman';
xlabel('Distance, km', 'FontSize', 25, 'FontName', 'Times New Roman');
ylabel('Depth, km', 'FontSize', 25, 'FontName', 'Times New Roman');
% axis equal;
% axis tight;
set(gca, 'YDir', 'reverse', 'FontSize', 25);
text(-0.12, 1.05, '(h)', 'FontSize', 30, 'FontName', 'Times New Roman', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'Units', 'normalized', 'Parent', gca, 'FontWeight', 'bold');
caxis([0 0.01]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot('Position', [0.7 0.15 0.25 0.8]);
[C, h] = contour(x, z, error_re, levels, 'LineWidth', 2);
c = colorbar;

caxis([min(levels) max(levels)]);
c.Label.String = 'Relative error, %';
c.Label.FontName='Times New Roman';
c.Label.FontSize=25;
c.FontSize = 25; 
c.FontName = 'Times New Roman';

colormap('jet');
xlabel('Distance, km', 'FontSize', 25, 'FontName', 'Times New Roman');
ylabel('Depth, km', 'FontSize', 25, 'FontName', 'Times New Roman');
text(-0.12, 1.05, '(i)', 'FontSize', 30, 'FontName', 'Times New Roman', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'Units', 'normalized', 'Parent', gca, 'FontWeight', 'bold');
clabel(C, h, 'FontSize', 10, 'Color', 'k', 'LabelSpacing', 300);
set(gca, 'YDir', 'reverse', 'FontSize', 25);
print("./tiff/homo_rsgfsm", '-dtiff', '-r500');