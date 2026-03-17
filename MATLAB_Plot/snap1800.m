nx0 = 539;
nz0 = 178;
nx = 501;
nz = 140;
dx = 10;
dz = 10;
x = [0:dx:(nx-1)*dx];
z = [0:dz:(nz-1)*dz];

% fid = fopen('/media/lqq/Data/wang/codes/Forward2D/output/snap800.bin', 'r');
fid = fopen('./snap/snap1800.bin', 'r'); %600*0.0005
snap0 = fread(fid, [nz0 nx0], 'float');
fclose(fid);
snap = snap0(20:159, 20:520); 
maxsnap = prctile(snap(:), 97);
minsnap = prctile(snap(:), 3);

fid = fopen('./mar/traveltime_FSM.dat', 'r');
time = fread(fid, [nz nx], 'float');
fclose(fid);
fid = fopen('./mar/traveltime_stagered_FSM.dat', 'r');
time1 = fread(fid, [nz nx], 'float');
fclose(fid);
time = time + 1.0/35;
time1 = time1 +1.0/35;

figure('Units', 'normalized', 'Position', [0 0 0.7 0.4]);
imagesc(x, z, snap);
colormap('gray');
% colorbar;
caxis([minsnap maxsnap]);
hold on;

contour(x,z,time, [0.9 0.9], 'g-','ShowText','off','LineWidth',2);
contour(x,z,time1, [0.9 0.9], 'r--','ShowText','off','LineWidth',2);
text(-0.08, 1, '(a)', 'FontSize', 20, 'FontName', 'Times New Roman', ...
     'Units', 'normalized', 'Parent', gca, 'FontWeight', 'bold');
xlabel('Distance(m)', 'FontSize', 20, 'FontName', 'Times New Roman');
ylabel('Depth(m)', 'FontSize', 20, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 20, 'Fontname', 'Times New Roman');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print('./picture/snap1800', '-dpng', '-r500');

