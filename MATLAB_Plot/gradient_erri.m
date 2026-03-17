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

fprintf("%f %f %f %f", min(ana(:)), max(ana(:)), min(fsm(:)), max(fsm(:)));
level = [0.1: 0.1: 0.9];

figure('units','normalized','position',[0,0,0.5,0.45]);

subplot('Position', [0.08 0.15 0.4 0.75]);
contour(x,z,ana,level,'k-','ShowText','off','LineWidth',1);
hold on;
contour(x,z,fsm,level,'g--','ShowText','off','LineWidth',1);
contour(x,z,fsm2nd,level,'b--','ShowText','off','LineWidth',1);
contour(x,z,sgfsm,level,'r--','ShowText','off','LineWidth',1);
hold off;

x_zoom = [50*dx, 80*dx, 80*dx, 50*dx, 50*dx];
z_zoom = [80*dz, 80*dz,50*dz, 50*dz, 80*dz];
hold on;
plot(x_zoom, z_zoom, 'b-', 'LineWidth', 1.5);
hold off;

xlabel('Distance, km', 'FontSize', 25, 'FontName', 'Times New Roman');
ylabel('Depth, km', 'FontSize', 25, 'FontName', 'Times New Roman');
text(-0.12, 1.1, '(a)', 'FontSize', 30, 'FontName', 'Times New Roman', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'Units', 'normalized', 'Parent', gca, 'FontWeight', 'bold');
set(gca, 'YDir', 'reverse', 'FontSize', 25);
legend('Analytical','FSM', 'Second-order FSM','RSGFSM', 'FontName', 'Times New Roman', ...
       'Location', 'south', 'Orientation', 'horizontal', 'FontSize', 13);

subplot('Position', [0.58 0.15 0.4 0.75]);
x_start_idx = 50;
x_end_idx = 80;
z_start_idx = 50;
z_end_idx = 80;

x_zoom = x(x_start_idx:x_end_idx);
z_zoom = z(z_start_idx:z_end_idx);
ana_zoom = ana(z_start_idx:z_end_idx, x_start_idx:x_end_idx);
fsm_zoom = fsm(z_start_idx:z_end_idx, x_start_idx:x_end_idx);
fsm2nd_zoom = fsm2nd(z_start_idx:z_end_idx, x_start_idx:x_end_idx);
sgfsm_zoom = sgfsm(z_start_idx:z_end_idx, x_start_idx:x_end_idx);

contour(x_zoom, z_zoom, ana_zoom, level, 'k-', 'ShowText','off','LineWidth',1.5);
hold on;
contour(x_zoom, z_zoom, fsm_zoom, level, 'g--', 'ShowText','off','LineWidth',1.5);
contour(x_zoom, z_zoom, fsm2nd_zoom, level, 'b--', 'ShowText','off','LineWidth',1.5);
contour(x_zoom, z_zoom, sgfsm_zoom, level, 'r--', 'ShowText','off','LineWidth',1.5);
hold off;

xlabel('Distance, km', 'FontSize', 25, 'FontName', 'Times New Roman');
ylabel('Depth, km', 'FontSize', 25, 'FontName', 'Times New Roman');
text(-0.12, 1.1, '(b)', 'FontSize', 30, 'FontName', 'Times New Roman', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
     'Units', 'normalized', 'Parent', gca, 'FontWeight', 'bold');
set(gca, 'YDir', 'reverse', 'FontSize', 25);
