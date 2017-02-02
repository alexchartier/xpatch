%% Matlab startup script
clear all

%% Set plotting defaults (edit as required)
S = get(0,'ScreenSize');
set(0,'DefaultFigurePosition',[8,82,S(3:4)/2])    % Figure position and size
set(0,'DefaultFigureColor',[1,1,1])               % Figure background color
set(0,'DefaultTextColor',[0,0,0])                 % Color for text and titles
set(0,'DefaultAxesColor',[1,1,1])              % Color fill for plot
set(0,'DefaultAxesYColor',[0,0,0])                % Color of y axis text and ticks
set(0,'DefaultAxesXColor',[0,0,0])                % Color of x axis text and ticks
set(0,'DefaultAxesZColor',[0,0,0])                % Color of z axis text and ticks
set(0,'DefaultTextFontSize',14)                    % Font size of titles
set(0,'DefaultAxesFontSize',14)                    % Font size of axes
set(0,'DefaultAxesFontName','arial')              % Axis Font
set(0,'DefaultTextFontName','arial')              % Axis Font
set(0,'defaultlinelinewidth',1)

%% Set path

% Sets the directory where matlab looks for startup.m
userpath(pwd)

% Plotting packages
path('../mexcdf/mexnc', path)
path('../mexcdf/snctools', path)


% Set paths to utilities here
UtilDir = '~/matlab_utils/';
path([UtilDir 'dart'], path)
path([UtilDir 'gps'], path)
path([UtilDir 'CalcMD5'], path)
path([UtilDir 'maths'], path)
path([UtilDir 'm_map'],path)
path([UtilDir 'geodetic'],path)
addpath([UtilDir 'AntarcticMappingTools'])
path([UtilDir 'export_fig'],path)
path([UtilDir 'midas'],path)
path([UtilDir 'raytrace'],path)

%% Environment variables (e.g. to external storage)
setenv('DataPath','/Volumes/Seagate/data/')

%% Set up fundamental constants
set_constants

