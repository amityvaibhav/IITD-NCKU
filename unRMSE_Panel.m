%%--------------ubRMSE Computation--------------%%

close all; clear; clc

%% =========================================================
%% Load datasets
%% =========================================================

load('WATERGAP_RZSM_1981-2020_India_Pentad_Percentile.mat','RZSM_pct')
WATERGAP_pct = RZSM_pct;

load('H08_RZSM_1981-2020_India_Pentad_Percentile.mat','RZSM_pct')
H08_pct = RZSM_pct;

load('ERA5L_RZSM_1981-2020_India_Pentad_Percentile.mat','RZSM_pct')
ERA5L_pct = RZSM_pct;

load('GLEAM_RZSM_1981-2020_India_Pentad_Percentile.mat','RZSM_pct')
GLEAM_pct = RZSM_pct;

clear RZSM_pct

%% =========================================================
%% Move TIME dimension first
%% =========================================================

moveTime = @(X) permute(X,[find(size(X)==max(size(X)),1) ...
                          setdiff(1:3,find(size(X)==max(size(X)),1))]);

WATERGAP_pct = moveTime(WATERGAP_pct);
H08_pct      = moveTime(H08_pct);
ERA5L_pct    = moveTime(ERA5L_pct);
GLEAM_pct    = moveTime(GLEAM_pct);

%% =========================================================
%% Harmonize spatial grid
%% =========================================================

lat = 129;
lon = 135;

WATERGAP_pct = WATERGAP_pct(:,1:lat,1:lon);
H08_pct      = H08_pct(:,1:lat,1:lon);
ERA5L_pct    = ERA5L_pct(:,1:lat,1:lon);
GLEAM_pct    = GLEAM_pct(:,1:lat,1:lon);

%% =========================================================
%% Fix latitude orientation
%% =========================================================

WATERGAP_pct = WATERGAP_pct(:,end:-1:1,:);
H08_pct      = H08_pct(:,end:-1:1,:);
ERA5L_pct    = ERA5L_pct(:,end:-1:1,:);
GLEAM_pct    = GLEAM_pct(:,end:-1:1,:);

%% =========================================================
%% Harmonize time dimension
%% =========================================================

T = 2919;

WATERGAP_pct = WATERGAP_pct(1:T,:,:);
H08_pct      = H08_pct(1:T,:,:);
ERA5L_pct    = ERA5L_pct(1:T,:,:);
GLEAM_pct    = GLEAM_pct(1:T,:,:);

fprintf('Final dataset size: %d time × %d lat × %d lon\n',T,lat,lon);

%% =========================================================
%% Reshape to time × grid
%% =========================================================

WATERGAP = reshape(WATERGAP_pct,T,[]);
H08      = reshape(H08_pct,T,[]);
ERA5L    = reshape(ERA5L_pct,T,[]);
GLEAM    = reshape(GLEAM_pct,T,[]);

N = size(WATERGAP,2);

%% =========================================================
%% Store datasets
%% =========================================================

Data  = {WATERGAP,H08,ERA5L,GLEAM};
Names = {'WATERGAP','H08','ERA5L','GLEAM'};

pairs = nchoosek(1:4,2);

%% =========================================================
%% Compute ubRMSE
%% =========================================================

UBRMSE = cell(6,1);

for p = 1:size(pairs,1)

    k = pairs(p,1);
    l = pairs(p,2);

    X = Data{k};
    Y = Data{l};

    % remove temporal mean
    Xc = X - mean(X,1,'omitnan');
    Yc = Y - mean(Y,1,'omitnan');

    % ubRMSE formula
    diff2 = (Xc - Yc).^2;
    ub = sqrt(mean(diff2,1,'omitnan'));

    ub_map = reshape(ub,lat,lon);

    UBRMSE{p} = ub_map;

    fname = [Names{k} '_vs_' Names{l} '_ubRMSE.mat'];
    save(fname,'ub_map','-v7.3')

    fprintf('Saved: %s\n',fname)

end

disp('ubRMSE maps computed successfully.')

%% =========================================================
%% Panel Plot (same style as Spearman figure)
%% =========================================================

%% India domain
lon0 = 65; lon1 = 96.5;
lat0 = 5;  lat1 = 36.5;

lon = linspace(lon0,lon1,lon);
lat = linspace(lat0,lat1,lat);

%% Grid settings
degStep   = 5;
gridColor = [0.5 0.5 0.5];
gridAlpha = 0.20;

%% Color limits for ubRMSE
CLIM = [0 25];

%% Figure
figure('Units','centimeters','Position',[2 2 32 18])

t = tiledlayout(2,3,'TileSpacing','compact','Padding','loose');
colormap(parula(256))

bgColor = [0.95 0.95 0.95];

applyGrid = @(ax) set(ax,...
    'XGrid','on','YGrid','on',...
    'GridColor',gridColor,'GridAlpha',gridAlpha,...
    'XTick',lon0:degStep:lon1,...
    'YTick',lat0:degStep:lat1,...
    'TickDir','out',...
    'Box','on',...
    'FontSize',11,...
    'FontWeight','bold');

titles = {'WATERGAP vs H08','WATERGAP vs ERA5 Land','WATERGAP vs GLEAM',...
          'H08 vs ERA5 Land','H08 vs GLEAM','ERA5 Land vs GLEAM'};

labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

%% Plot panels
for i = 1:6

    nexttile

    img = imagesc(lon,lat,UBRMSE{i});
    set(img,'AlphaData',~isnan(UBRMSE{i}));

    set(gca,'YDir','normal','Color',bgColor)

    axis image
    axis tight

    clim(CLIM)

    set(gca,'LineWidth',1.2)

    ax = gca; applyGrid(ax);

    title(titles{i},'FontWeight','bold')

    text(0.02,0.98,labels{i},'Units','normalized',...
        'FontSize',12,'FontWeight','bold','VerticalAlignment','top');

end

%% Vertical colorbar
drawnow

ax = findall(gcf,'Type','axes');
ax = flipud(ax);

axTop = ax(3);
axBot = ax(6);

posTop = axTop.Position;
posBot = axBot.Position;

bottom = posBot(2);
top    = posTop(2) + posTop(4);

height = top-bottom;

left   = posTop(1) + posTop(3) + 0.02;
width  = 0.018;

cb = colorbar('Position',[left bottom width height]);

cb.Ticks = 0:5:25;
cb.FontSize = 13;
cb.FontWeight = 'bold';

cb.Label.String = 'Unbiased RMSE (percentile)';
cb.Label.FontWeight = 'bold';
cb.Label.FontSize = 14;

%% Export
set(gcf,'Color','w')
drawnow

exportgraphics(gcf,'output_file.tif','Resolution',600);

disp('ubRMSE panel exported successfully.')
