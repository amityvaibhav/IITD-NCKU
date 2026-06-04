%%====================================================================%%
% FD_India: Structural Spread (IQR) Computation and Visualization
%
% Computes structural spread of flash drought frequency across
% WATERGAP, H08, ERA5-Land, and GLEAM
%
% Author  : Prajwal Giri
% Project : FD_India
%%====================================================================%%

close all; clear; clc

fprintf('\nFD_India – Structural Spread Computation\n');
fprintf('Run date: %s\n\n',datestr(now));

%% =========================================================
%% Load flash drought frequency datasets
%% =========================================================

load('WATERGAP_FD_frequency_total_India.mat','FD_freq_total'); WG = FD_freq_total;
load('H08_FD_frequency_total_India.mat','FD_freq_total');      H08 = FD_freq_total;
load('ERA5L_FD_frequency_total_India.mat','FD_freq_total');    ERA = FD_freq_total;
load('GLEAM_FD_frequency_total_India.mat','FD_freq_total');    GLEAM = FD_freq_total;

clear FD_freq_total

fprintf('Datasets loaded successfully\n');

%% =========================================================
%% Harmonize spatial grid (important for GLEAM mismatch)
%% =========================================================

lat = min([size(WG,1), size(H08,1), size(ERA,1), size(GLEAM,1)]);
lon = min([size(WG,2), size(H08,2), size(ERA,2), size(GLEAM,2)]);

WG    = WG(1:lat,1:lon);
H08   = H08(1:lat,1:lon);
ERA   = ERA(1:lat,1:lon);
GLEAM = GLEAM(1:lat,1:lon);

fprintf('Spatial grid harmonized: %d lat × %d lon\n',lat,lon);

%% =========================================================
%% FIX 180° FLIP (latitude orientation)
%% =========================================================

WG    = WG(end:-1:1,:);
H08   = H08(end:-1:1,:);
ERA   = ERA(end:-1:1,:);
GLEAM = GLEAM(end:-1:1,:);

fprintf('Latitude orientation corrected\n');

%% =========================================================
%% Stack products
%% dimension: lat × lon × product
%% =========================================================

DATA = cat(3,WG,H08,ERA,GLEAM);

%% =========================================================
%% Compute quartiles
%% =========================================================

Q25 = prctile(DATA,25,3);
Q75 = prctile(DATA,75,3);

%% =========================================================
%% Structural spread (IQR)
%% =========================================================

IQR_map = Q75 - Q25;

fprintf('Structural spread (IQR) computed\n');

%% =========================================================
%% Save spread dataset
%% =========================================================

save('FD_Frequency_StructuralSpread_India.mat',...
     'IQR_map','Q25','Q75','-v7.3');

fprintf('Spread dataset saved\n');

%% =========================================================
%% Visualization
%% =========================================================

lon0 = 65; lon1 = 96.5;
lat0 = 5;  lat1 = 36.5;

lonVec = linspace(lon0,lon1,lon);
latVec = linspace(lat0,lat1,lat);

%% Mask tiny values
maskThr = 0.001;
IQR_map(IQR_map < maskThr) = NaN;

%% Color scaling
vals = IQR_map(~isnan(IQR_map));
CLIM = [0 prctile(vals,98)];

%% Plot map
figure('Units','centimeters','Position',[2 2 15 10])
set(gcf,'Color','w')

colormap(turbo(256))

img = imagesc(lonVec,latVec,IQR_map);
set(img,'AlphaData',~isnan(IQR_map));

set(gca,'YDir','normal','Color','w')

axis image
axis tight
clim(CLIM)

set(gca,'LineWidth',1.2)

%% Grid styling
degStep   = 5;
gridColor = [0.6 0.6 0.6];
gridAlpha = 0.15;

set(gca,...
    'XGrid','on','YGrid','on',...
    'GridColor',gridColor,...
    'GridAlpha',gridAlpha,...
    'XTick',lon0:degStep:lon1,...
    'YTick',lat0:degStep:lat1,...
    'TickDir','out',...
    'FontSize',11,...
    'FontWeight','bold');

title('Structural Spread (IQR of Flash Drought Frequency)','FontWeight','bold')

%% Colorbar
cb = colorbar;
cb.FontSize = 13;
cb.FontWeight = 'bold';
cb.Label.String = 'Structural Spread (IQR)';
cb.Label.FontWeight = 'bold';

%% Export figure
exportgraphics(gcf,...
'RZSM_StructuralSpread_India_final.tif',...
'Resolution',600);

fprintf('Structural spread figure exported successfully\n');