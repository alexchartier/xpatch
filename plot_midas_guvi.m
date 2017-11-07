%% plot_midas_guvi.m

%% Set input parameters
IPath = '~/xpatch/data/midas_output/outPC10min3{yymmmdd-HHMM}.mat';

Time = datenum(2015, 12, 20, 18, 20, 0);
before = 20/60/24;
after = 20/60/24;
lat_cutoff = 45; 
crd = 'geo';
fname = '~/Downloads/GUVI_sp_v010r00_2015354_REV76053.L1B';

%% load MIDAS
D = tec(load(filename(IPath, Time)));

if strcmp(crd, 'mag')
   Sph = cartsph([D.X(:), D.Y(:), D.Z(:)] * geomag);
else
   Sph = cartsph([D.X(:)'; D.Y(:)'; D.Z(:)']');
end
Lat = reshape(rad2deg(Sph(:, 2)), [length(D.Lat), length(D.Lon)]);
Lon = reshape(rad2deg(Sph(:, 3)), [length(D.Lat), length(D.Lon)]);

%% Load GUVI
t_ms = double(ncread(fname, 'TIME'));
wavelengths = double(ncread(fname, 'WAVELENGTHS'));
pixel_lat = double(ncread(fname, 'PIXELLATITUDE'));
pixel_lon = double(ncread(fname, 'PIXELLONGITUDE'));
spectra = double(ncread(fname, 'PIXELSPECTRA'));

t = datenum(2015, 1, 354, 0, 0, t_ms / 1000);
ind_1356 = repmat(wavelengths >= 1.355E3 & wavelengths <= 1.357E3, [1, 1, 1937]);


%% Plot MIDAS
clf
subplot(2, 1, 1)
colormap jet
hold on
m_proj('Stereographic', 'latitude', 90, 'radius', 40, 'rotation', 315)

[~, hC] = m_contourf(Lon, Lat, squeeze(D.F), 100);
set(hC,'LineStyle','none')
m_coast('color', 'w', 'LineWidth', 2)
title(filename('{dd/mm/yyyy HH:MM UT}                                            ', D.Time))
caxis([0 20]) % side colour axis
h = colorbar;
ylabel(h, 'TEC (TECU)')

m_grid('color', 'k', 'FontSize', 20)  % TODO: ticks every 10 degrees
set(gca, 'FontSize', 20)
title(filename('{yyyy/mm/dd HH:MM}UT', Time))

% Overplot GUVI loc
pixel_lon(pixel_lon > 180) = pixel_lon(pixel_lon > 180) - 360;
blackline = m_plot(pixel_lon(:), pixel_lat(:), 'kx');

%% Plot GUVI data
subplot(2, 1, 2)
t_2d = repmat(t, [1, 14])';
tind = t_2d > (Time - before) & t_2d < (Time + after);
scatter(pixel_lon(tind), pixel_lat(tind), 50, spectra(tind), 'filled'); 
xlabel('Lon')
ylabel('Lat')
set(gca, 'clim', [0, 100])
% 
% %% plot GUVI locs
% t_2d = repmat(t, [1, 14])';
% datestr(t_2d(tind))
