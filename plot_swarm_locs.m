%% plot_swarm_locs.m
% Script to plot Swarm satellite locations around the time of the eclipse

%% Set input parameters
RootDir = '~/xpatch/data/eclipse/'; 
Sats = ['A', 'B'];

isr_locs = [67, -50; 65.1, -147.4; 44.15, -72.4; 18.4, -66.6; -11.9, -76.9];
Time = datenum(2017, 8, 25);
m_proj('miller');
m_coast('color', 'k', 'linewidth', 2);
m_grid('linestyle','none','tickdir','out');
colors = ['g', 'm'];
hold on
ct = 0;
for Sat = Sats
    ct = ct + 1;
    SwarmPath = [RootDir, 'lp/SW_OPER_EFI', Sat, ...
        '_PL_1B_{yyyymmdd}T000000_{yyyymmdd}T235959_0403.CDF/SW_OPER_EFI', Sat, ...
        '_PL_1B_{yyyymmdd}T000000_{yyyymmdd}T235959_0403_MDR_EFI_PL.cdf'];
    %% load
    lats = cell2mat(cdfread(filename(SwarmPath, Time), 'Variable', 'Latitude'));
    lons = cell2mat(cdfread(filename(SwarmPath, Time), 'Variable', 'Longitude'));
    
    
    m_plot(lons, lats, ['.', colors(ct)])
end
h = m_plot(isr_locs(:, 2), isr_locs(:, 1), 'pr', 'markersize', 40, ...
    'markerfacecolor', 'r');

h = m_plot(-72.4, 44.15, 'pr', 'markersize', 40, ...
    'markerfacecolor', 'b');

uistack(h, 'top')

hold off











