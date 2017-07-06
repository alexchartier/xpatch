%% plot_midas_tec.m
% Script to plot the TEC output of MIDAS from a polar perspective

%% Set input parameters
RootDir = './'; 
IPath = [RootDir, 'data/test/outPC10min3{yymmmdd-HHMM}.mat'];
Sat = 'B';
SwarmPath = [RootDir, 'lp/SW_OPER_EFI', Sat, ...
   '_PL_1B_{yyyymmdd}T000000_{yyyymmdd}T235959_0403_MDR_EFI_PL.cdf'];

Time = datenum(2016, 10, 12, 12, 50, 0);
before = 20/60/24;
after = 10/60/24;
lat_cutoff = 50;
% crd = 'mag';
crd = 'geo';

%% load
D = tec(load(filename(IPath, Time)));

if strcmp(crd, 'mag')
   Sph = cartsph([D.X(:), D.Y(:), D.Z(:)] * geomag);
else
   Sph = cartsph([D.X(:)'; D.Y(:)'; D.Z(:)']');
end
Lat = reshape(rad2deg(Sph(:, 2)), [length(D.Lat), length(D.Lon)]);
Lon = reshape(rad2deg(Sph(:, 3)), [length(D.Lat), length(D.Lon)]);

%% Plot
colormap jet
hold on
m_proj('Stereographic', 'latitude', 90, 'radius', 40, 'rotation', 90)

[~, hC] = m_contourf(Lon, Lat, squeeze(D.F), 100);
set(hC,'LineStyle','none')
m_coast('color', 'w', 'LineWidth', 2)
title(filename('{dd/mm/yyyy HH:MM UT}                                            ', D.Time))
caxis([0 15]) % side colour axis
h = colorbar;
ylabel(h, 'TEC (TECU)')

m_grid('color', 'k', 'FontSize', 20)  % TODO: ticks every 10 degrees
set(gca, 'FontSize', 20)
title(filename('{yyyy/mm/dd HH:MM}UT', Time))
% % blackline = m_plot(Swarm_lon, Swarm_lat, 'kx');
% blackline = m_plot([305 323], [62, 51], 'kx-');
% redline = m_plot([322 347], [61, 50], 'rx-');
% 
% uistack(blackline, 'top')
hold off
















