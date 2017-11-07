%% plot_rbsp_tec.m
% Script to plot the TEC output of MIDAS from a polar perspective

%% Set input parameters
RootDir = '/Volumes/Seagate/data/swarm/';
IPath = [RootDir, 'midas/outPC10min3{yymmmdd-HHMM}.mat'];
RBSPPath = '~/xpatch/Oxygen_20150623_info.dat';
OPath = '~/rbsp/{yyyymmmdd-HHMM}.png';
Times = datenum(2015, 6, 23, 3, 0, 0):10/60/24:datenum(2015, 6, 23, 6, 0, 0);
latlim = [55, 65];
lonlim = [300, 330];

crd = 'mag';
t = 1;

%% Loop over times
for t = 1:length(Times)
   %% load TEC 
   D = tec(load(filename(IPath, Times(t))));
   
   if strcmp(crd, 'mag')
      Sph = cartsph([D.X(:), D.Y(:), D.Z(:)] * geomag);
   else
      Sph = cartsph([D.X(:)'; D.Y(:)'; D.Z(:)']');
   end
   Lat = reshape(rad2deg(Sph(:, 2)), [length(D.Lat), length(D.Lon)]);
   Lon = reshape(rad2deg(Sph(:, 3)), [length(D.Lat), length(D.Lon)]);
   Lon(Lon < 0) = Lon(Lon < 0) + 360;
   
   %% Load RBSP position
   Txt = asciiread(RBSPPath);
   Vals = str2num(Txt(4:end, :));
   RBSP.Time = datenum(Vals(:, 1), Vals(:, 2), Vals(:, 3), Vals(:, 4), Vals(:, 5), 0);
   RBSP.Mlat = Vals(:, 8);
   RBSP.Mlon = Vals(:, 10);
   
   %% Plot
   clf
   colormap jet
   hold on
   m_proj('Stereographic', 'latitude', 90, 'radius', 40, 'rotation', 270)
   
   [~, hC] = m_contourf(Lon, Lat, squeeze(D.F), 100);
   set(hC,'LineStyle','none')
   % m_coast('color', 'w', 'LineWidth', 1.2)
   % title(filename('{dd/mm/yyyy HH:MM UT}                                            ', D.Time))
   caxis([0 20]) % side colour axis
   h = colorbar;
   ylabel(h, 'TEC (TECU)')
   
   m_grid('color', 'k', 'FontSize', 20)
   set(gca, 'FontSize', 20)
   title(filename('{yyyy/mm/dd HH:MM}UT', Times(t)))
   foot_mlon = RBSP.Mlon(RBSP.Time == Times(t));
   foot_mlat = RBSP.Mlat(RBSP.Time == Times(t));
   fprintf('%2.2f %2.2f', foot_mlon, foot_mlat)
   redline = m_plot(foot_mlon, foot_mlat, 'rx', 'markersize', 25, 'linewidth', 5);
   
   uistack(redline, 'top')
   hold off
   
   %%
   export_fig(filename(OPath, Times(t)))
   close
end
















