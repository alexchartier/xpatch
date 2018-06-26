%% plot_midas_tec.m
% Script to plot the TEC output of MIDAS from a polar perspective

%% Set input parameters
% RootDir = '~/xpatch/data/';
clear
RootDir = '~/midas-3/matlab/scripts/data/ningchao/';
IPath = [RootDir, '/midas_{yymmmdd-HHMM}.mat'];
Sat = 'B';
SwarmPath = ['/Volumes/Seagate/data/swarm/lp/SW_EXTD_EFI', Sat, ...
   '_LP_HM_{yyyymmdd}T000000_{yyyymmdd}T235959_0101.cdf'];

Time = datenum(2015, 3, 17, 18, 0, 0);
before = 15/60/24;
after = 25/60/24;
lat_cutoff = 45; 
crd = 'geo';
% crd = 'mag';


%% load
D = load(filename(IPath, Time));
% D = tec(D1);

if strcmp(crd, 'mag')
   Sph = cartsph([D.X(:), D.Y(:), D.Z(:)] * geomag);
else
   Sph = cartsph([D.X(:)'; D.Y(:)'; D.Z(:)']');
end

% Shape = [length(D.Lat), length(D.Lon)];
Shape = [size(D.Lat, 2), size(D.Lon, 3)];
Lat = reshape(rad2deg(Sph(:, 2)), Shape);
Lon = reshape(rad2deg(Sph(:, 3)), Shape);


%% Load Swarm
Swarm.Lat = cell2mat(cdfread(filename(SwarmPath, Time), 'Variable', {'Latitude'}));
Swarm.Lon = cell2mat(cdfread(filename(SwarmPath, Time), 'Variable', {'Longitude'}));
Swarm.Timestamp = cdfread(filename(SwarmPath, Time), 'Variable', {'Timestamp'});
Swarm.Ne = cdfread(filename(SwarmPath, Time), 'Variable', {'n'});

i = zeros(size(Swarm.Timestamp));
for i = 1:length(Swarm.Timestamp)
    Swarm.Time(i) = todatenum(Swarm.Timestamp{i});
end

tind = Swarm.Time > Time - before & Swarm.Time < Time + after;
ind = tind & Swarm.Lat' > 50;
Swarm_lon = Swarm.Lon(ind);
Swarm_lat = Swarm.Lat(ind);
Swarm_t = Swarm.Time(ind);
Swarm_Ne = cell2mat(Swarm.Ne(ind));

%% Plot
figure
%colormap jet
hold on
m_proj('Stereographic', 'latitude', 90, 'radius', 40, 'rotation', 222)

[~, hC] = m_contourf(Lon, Lat, squeeze(D.F), 100);
set(hC,'LineStyle','none')
m_coast('color', 'w', 'LineWidth', 2)
title(filename('{dd/mm/yyyy HH:MM UT}                                            ', D.Time))
caxis([0 16]) % side colour axis
h = colorbar;
ylabel(h, 'TEC (TECU)')

m_grid('color', 'k', 'FontSize', 20)  % TODO: ticks every 10 degrees
set(gca, 'FontSize', 20)
title(['Swarm ', Sat, ...
    filename(' {yyyy/mm/dd HH:MM}UT', Time)])
blackline = m_plot(Swarm_lon, Swarm_lat, 'kx');

mins = [ceil(min(Swarm_t * 1440)):2:floor(max(Swarm_t * 1440))] / 1440;
Swarm_mins = zeros(size(mins));
for t = 1:length(mins)
    Swarm_mins(t) = closest(Swarm_t, mins(t));
end
minind = ismember(Swarm_t, Swarm_mins);
redline = m_plot(Swarm_lon(minind), Swarm_lat(minind), 'r.', 'markersize', 20);

hold off


% saveas(gcf, filename('~/Dropbox/papers/patch_grl/midas_plots/{yyyymmmdd-HHMM}.jpg', Time))
%%
figure; 
plot(Swarm_t, Swarm_Ne, 'LineWidth', 3) 
axis tight
ylim([0 max(Swarm_Ne) + 5E4])
xticks(Swarm_mins)
datetick('keepticks')
grid on
grid minor
xlabel('Time (UT)')
ylabel('Electrons / cm^3')
set(gca, 'LineWidth', 2, 'FontSize', 23)













