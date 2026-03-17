nx = 101;
nz = 101;
dx = 0.01;
dz = 0.01;
x = [0:dx:(nx-1)*dx];
z = [0:dz:(nz-1)*dz];

fid=fopen(['./gradient/traveltime_gradient_analy.bin'],'r');
ana=fread(fid,[nz nx],'float');
fclose(fid);

fid=fopen(['./gradient/traveltime_FSM.dat'],'r');
fsm=fread(fid,[nz nx],'float');
fclose(fid);

fid=fopen(['./gradient/traveltime_2ndFSM.dat'],'r');
fsm2nd=fread(fid,[nz nx],'float');
fclose(fid);

fid=fopen(['./gradient/traveltime_stagered_FSM.dat'],'r');
sgfsm=fread(fid,[nz nx],'float');
fclose(fid);

error_abs1 = abs(fsm-ana);
error_abs2 = abs(fsm2nd-ana);
error_abs3 = abs(sgfsm-ana);
fprintf("%f %f %f", max(error_abs1(:)), max(error_abs2(:)), max(error_abs3(:)));


figure('units','normalized','position',[0,0,0.9,0.45]);

subplot('Position', [0.05 0.12 0.24 0.8]);
imagesc(x, z, error_abs1);
colormap('jet');
c = colorbar;
c.Label.String = 'Absolute error, s';
c.Label.FontName='Times New Roman';
c.Label.FontSize=25;
c.FontSize = 25; 
c.FontName = 'Times New Roman';
xlabel('Distance, km', 'FontSize', 25, 'FontName', 'Times New Roman');
ylabel('Depth, km', 'FontSize', 25, 'FontName', 'Times New Roman');
axis equal;
axis tight;
text(-0.1, 1.1, '(a)', 'FontSize', 30, 'FontName', 'Times New Roman', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'Units', 'normalized', 'Parent', gca, 'FontWeight', 'bold');
set(gca, 'FontSize', 25);
caxis([0 0.01]);

subplot('Position', [0.375 0.12 0.24 0.8]);
imagesc(x, z, error_abs2);
colormap('jet');
c = colorbar;
c.Label.String = 'Absolute error, s';
c.Label.FontName='Times New Roman';
c.Label.FontSize=25;
c.FontSize = 25; 
c.FontName = 'Times New Roman';
xlabel('Distance, km', 'FontSize', 25, 'FontName', 'Times New Roman');
ylabel('Depth, km', 'FontSize', 25, 'FontName', 'Times New Roman');
axis equal;
axis tight;
text(-0.1, 1.1, '(b)', 'FontSize', 30, 'FontName', 'Times New Roman', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'Units', 'normalized', 'Parent', gca, 'FontWeight', 'bold');
set(gca, 'FontSize', 25);
caxis([0 0.01]);

subplot('Position', [0.7 0.12 0.24 0.8]);
imagesc(x, z, error_abs3);
colormap('jet');
c = colorbar;
c.Label.String = 'Absolute error, s';
c.Label.FontName='Times New Roman';
c.Label.FontSize=25;
c.FontSize = 25; 
c.FontName = 'Times New Roman';
xlabel('Distance, km', 'FontSize', 25, 'FontName', 'Times New Roman');
ylabel('Depth, km', 'FontSize', 25, 'FontName', 'Times New Roman');
axis equal;
axis tight;
text(-0.1, 1.1, '(c)', 'FontSize', 30, 'FontName', 'Times New Roman', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'Units', 'normalized', 'Parent', gca, 'FontWeight', 'bold');
set(gca, 'FontSize', 25);
caxis([0 0.01]);

% print("./picture/gradient_abs_error", '-dpng', '-r500');