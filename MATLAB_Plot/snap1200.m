nx0 = 539;
nz0 = 178;
nx = 501;
nz = 140;
dx = 0.01;
dz = 0.01;
x = [0:dx:(nx-1)*dx];
z = [0:dz:(nz-1)*dz];

% fid = fopen('/media/lqq/Data/wang/codes/Forward2D/output/snap800.bin', 'r');
fid = fopen('./snap/snap1200.bin', 'r'); %600*0.0005
snap0 = fread(fid, [nz0 nx0], 'float');
fclose(fid);
snap = snap0(20:159, 20:520); 
maxsnap = prctile(snap(:), 99);
minsnap = prctile(snap(:), 1);

fid = fopen('./mar/traveltime_FSM.dat', 'r');
time = fread(fid, [nz nx], 'float');
fclose(fid);
fid = fopen('./mar/traveltime_stagered_FSM.dat', 'r');
time1 = fread(fid, [nz nx], 'float');
fclose(fid);
time = time + 1.0/35;
time1 = time1 +1.0/35;

figure('units','normalized','position',[0.05,0.05,0.4,0.3]);
imagesc(x, z, snap);
colormap('gray');
% colorbar;
caxis([minsnap maxsnap]);
hold on;

contour(x,z,time, [0.6 0.6], 'g-','ShowText','off','LineWidth',3);
contour(x,z,time1, [0.6 0.6], 'r--','ShowText','off','LineWidth',3);
text(-0.1, 1, '(a)', 'FontSize', 30, 'FontName', 'Times New Roman', ...
     'Units', 'normalized', 'Parent', gca, 'FontWeight', 'bold');
xlabel('Distance, km', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Depth, km', 'FontSize', 20, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 20, 'Fontname', 'Times New Roman');
set(gca, 'YDir', 'reverse', 'FontSize', 25);
legend('FSM', 'RSGFSM', 'Location', 'northeast', 'FontName', 'Times New Roman', 'FontSize', 25);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print('./picture/snap1200', '-dpng', '-r500');

