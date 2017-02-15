%% plot_midas_tec.m
% Script to plot the TEC output of MIDAS from a polar perspective

%% Set input parameters
RootDir = '/Volumes/Seagate/data/swarm/'; 
IPath = [RootDir, 'midas/outPC10min3{yymmmdd-HHMM}.mat'];
Sat = 'B';
SwarmPath = [RootDir, 'lp/SW_OPER_EFI', Sat, ...
   '_PL_1B_{yyyymmdd}T000000_{yyyymmdd}T235959_0403_MDR_EFI_PL.cdf'];

Time = datenum(2015, 12, 20, 16, 40, 0);
Window = 15/60/24;
lat_cutoff = 55; 

%% load
D = tec(load(filename(IPath, Time)));
Sph = cartsph([D.X(:)'; D.Y(:)'; D.Z(:)']');
Lat = reshape(rad2deg(Sph(:, 2)), [length(D.Lat), length(D.Lon)]);
Lon = reshape(rad2deg(Sph(:, 3)), [length(D.Lat), length(D.Lon)]);

Swarm_lat = cell2mat(cdfread(filename(SwarmPath, Time), 'variables', 'Latitude'));
Swarm_lon = cell2mat(cdfread(filename(SwarmPath, Time), 'variables', 'Longitude'));
Swarm_ne = cell2mat(cdfread(filename(SwarmPath, Time), 'variables', 'n'));
Swarm_ne_err = cell2mat(cdfread(filename(SwarmPath, Time), 'variables', 'n_error'));
Swarm_time = cdfread(filename(SwarmPath, Time), 'variables', 'Timestamp');

%% Swarm preprocessing
Swarm_t = nan(length(Swarm_time), 1);
for i = 1:length(Swarm_time)
   Swarm_t(i) = todatenum(Swarm_time{i});
end

latind = Swarm_lat > lat_cutoff;
tind = Swarm_t > Time - Window & Swarm_t < Time + Window;

Swarm_t = Swarm_t(latind & tind);
Swarm_lon = Swarm_lon(latind & tind);
Swarm_lat = Swarm_lat(latind & tind);
Swarm_ne = Swarm_ne(latind & tind);
Swarm_ne_err = Swarm_ne_err(latind & tind);

%% Plot
subplot(2, 1, 1)
colormap jet
hold on
m_proj('Stereographic', 'latitude', 90, 'radius', 35, 'rotation', 200)

[~, hC] = m_contourf(Lon, Lat, squeeze(D.F), 100);
set(hC,'LineStyle','none')
m_coast('color', 'w', 'LineWidth', 1.2)
title(filename('{dd/mm/yyyy HH:MM UT}                                            ', D.Time))
caxis([0 14]) % side colour axis
colorbar

m_grid('color', 'k', 'FontSize', 20)
set(gca, 'FontSize', 20)
blackline = m_plot(Swarm_lon, Swarm_lat, 'kx');
uistack(blackline, 'top')
hold off

subplot(2, 1, 2)
plot(Swarm_t, Swarm_ne)
xlim([min(Swarm_t), max(Swarm_t)])
ylabel('Electron density (cm^-^3)')
datetick('keeplimits')
grid on
grid minor
















