%% plot_midas_tec.m
% Script to plot the TEC output of MIDAS from a polar perspective

%% Set input parameters
RootDir = '/Volumes/Seagate/data/swarm/'; 
IPath = [RootDir, 'midas/outPC10min3{yymmmdd-HHMM}.mat'];
Sat = 'B';
SwarmPath = [RootDir, 'lp/SW_OPER_EFI', Sat, ...
   '_PL_1B_{yyyymmdd}T000000_{yyyymmdd}T235959_0403_MDR_EFI_PL.cdf'];

Time = datenum(2015, 12, 20, 17, 0, 0);
lat_cutoff = 55; 


%% load
D = tec(load(filename(IPath, Time)));
Sph = cartsph([D.X(:)'; D.Y(:)'; D.Z(:)']');
Lat = reshape(rad2deg(Sph(:, 2)), [length(D.Lat), length(D.Lon)]);
Lon = reshape(rad2deg(Sph(:, 3)), [length(D.Lat), length(D.Lon)]);

Swarm_lat = cell2mat(cdfread(filename(SwarmPath, Time), 'variables', 'Latitude'));
Swarm_lon = cell2mat(cdfread(filename(SwarmPath, Time), 'variables', 'Longitude'));

Swarm_time = cdfread(filename(SwarmPath, Time), 'variables', 'Timestamp');
Swarm_t = nan(length(Swarm_time));
for i = 1:length(Swarm_time)
   Swarm_t(i) = todatenum(Swarm_time{i});
end

latind = Swarm_lat > lat_cutoff;
Swarm_t = Swarm_t(latind);
Swarm_lon = Swarm_lon(latind);
Swarm_lat = Swarm_lat(latind);


%% Plot
m_proj('Stereographic', 'latitude', 90, 'radius', 35, 'rotation', -30)
m_pcolor(Lon, Lat, squeeze(D.F));
% set(hC,'LineStyle','none')
title(datestr(D.Time))
caxis([0 12]) % side colour axis
colorbar

for i = 1:length(Swarm_lat)
  m_plot(Swarm_lon(i), Swarm_lat(i), 'x')
end

m_coast('color', 'w')
m_grid('color', 'r', 'FontSize', 24)
set(gca, 'FontSize', 24)
