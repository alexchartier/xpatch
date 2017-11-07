%% plot_swarm_locs.m
% Script to plot Swarm satellite locations around the time of the eclipse

%% Set input parameters
RootDir = '~/xpatch/data/eclipse/'; 
Sat = 'A';
% SwarmPath = [RootDir, 'gps/SW_OPER_TEC', Sat, ...
%    'TMS_2F_{yyyymmdd}T000000_{yyyymmdd}T235959_0202.cdf'];

SwarmPath = [RootDir, 'lp/SW_OPER_EFI', Sat, ...
   '_PL_1B_{yyyymmdd}T000000_{yyyymmdd}T235959_0403.CDF/SW_OPER_EFI', Sat, ...
   '_PL_1B_{yyyymmdd}T000000_{yyyymmdd}T235959_0403_MDR_EFI_PL.cdf'];

lat_bounds = [15, 70];
lon_bounds = [-135, -65];

%% Time loop 
clf
m_proj('miller','lat',[0 82], 'lon', [-180 -50]);
m_coast('color',[0 .6 0]);
hold on
Time = datenum(2017, 8, 18, 17, 0, 0);
plus = 3/24;
endtime = datenum(2017, 8, 25);
colours = {'r', 'g', 'b', 'k','c', 'm',  'y'};
timelist = [];
handles = [];
count = 1;
while Time < endtime
    timelist = [timelist, Time];
    %% load
    lats = cell2mat(cdfread(filename(SwarmPath, Time), 'Variable', 'Latitude'));
    lons = cell2mat(cdfread(filename(SwarmPath, Time), 'Variable', 'Longitude'));
    times_obj = cdfread(filename(SwarmPath, Time), 'Variable', 'Timestamp');
    times = zeros(length(times_obj), 1);
    for t = 1:length(times_obj)
        times(t) = todatenum(times_obj{t});
    end
    
    
    %% Indices
    timeind = (times > Time) & (times < Time + plus);
    latind = (lats > lat_bounds(1)) & (lats < lat_bounds(2));
    lonind = (lons > lon_bounds(1)) & (lons < lat_bounds(2));
    ind = timeind & latind & lonind;
    
    
    %% Plot
    handles = [handles, m_plot(lons(ind), lats(ind), ['.', colours{count}])];
    
    
    Time = Time + 1;
    count = count + 1;
end

timestrs = cellstr(datestr(timelist, 'mm-dd'));
m_legend(handles, datestr(timestrs{1}),datestr(timestrs{2}),datestr(timestrs{3}),...
        datestr(timestrs{4}),datestr(timestrs{5}),datestr(timestrs{6}),...
        datestr(timestrs{7}));
    
    %%
m_grid('linestyle','none','box','fancy','tickdir','out');











