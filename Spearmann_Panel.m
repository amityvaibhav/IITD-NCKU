%%-------------- Spearman Correlation Panel--------------%%
close all; clear; clc

%%---------------- Provide input files -----------------------------%%
load('file01.mat','rho_map');        WG_H08 = rho_map;

%%---------------- India domain (0.25°) ----------------------------%%
lon0 = 65; lon1 = 96.5;
lat0 = 5;  lat1 = 36.5;

[nr,nc] = size(WG_H08);

lon = linspace(lon0, lon1, nc);
lat = linspace(lat0, lat1, nr);

%%---------------- Grid settings -----------------------------------%%
degStep   = 5;
gridColor = [0.5 0.5 0.5];
gridAlpha = 0.20;

%%---------------- Color limits ------------------------------------%%
CLIM = [0 1];

%%---------------- Figure layout -----------------------------------%%
figure('Units','centimeters','Position',[2 2 32 18]);

t = tiledlayout(2,3,'TileSpacing','compact','Padding','loose');
colormap((parula(256)));

bgColor = [0.95 0.95 0.95];

applyGrid = @(ax) set(ax,...
    'XGrid','on','YGrid','on',...
    'GridColor',gridColor,'GridAlpha',gridAlpha,...
    'XTick',lon0:degStep:lon1,...
    'YTick',lat0:degStep:lat1,...
    'TickDir','out',...
    'Box','on',...
    'FontSize',10,...
    'FontWeight','bold');

DATA = {WG_H08, WG_ERA, WG_GLEAM, H08_ERA, H08_GLEAM, ERA_GLEAM};
titles = {'WATERGAP vs H08','WATERGAP vs ERA5 Land','WATERGAP vs GLEAM',...
          'H08 vs ERA5 Land','H08 vs GLEAM','ERA5 Land vs GLEAM'};

labels = {'(a)','(b)','(c)','(d)','(e)','(f)'};

%%---------------- Plot panels -------------------------------------%%
for i = 1:6
    
    nexttile
    
    img = imagesc(lon,lat,DATA{i});
    set(img,'AlphaData',~isnan(DATA{i}));
    
    set(gca,'YDir','normal','Color',bgColor)
    axis tight
    axis image
    
    clim(CLIM)
    
    set(gca,'LineWidth',1.2)
    
    ax = gca; applyGrid(ax);
    
    title(titles{i})
    
    text(0.02,0.98,labels{i},'Units','normalized',...
        'FontSize',12,'FontWeight','bold','VerticalAlignment','top');

end

%%---------------- Shared vertical colorbar ------------------------%%
drawnow

ax = findall(gcf,'Type','axes');
ax = flipud(ax);

% use right column axes for vertical span
ax3 = ax(3);
ax6 = ax(6);

pos3 = ax3.Position;
pos6 = ax6.Position;

bottom = pos6(2);
top    = pos3(2) + pos3(4);

height = top - bottom;

left   = pos3(1) + pos3(3) + 0.02;
width  = 0.018;

cb = colorbar('Position',[left bottom width height]);
cb.Orientation = 'vertical';

% Spearman scale
cb.Ticks = 0:0.2:1;

% Make colorbar numbers bigger and bold
cb.FontSize = 13;
cb.FontWeight = 'bold';

cb.Label.String = 'Spearman Rank Correlation (\rho)';
cb.Label.FontWeight = 'bold';
cb.Label.FontSize = 14;

%%---------------- Export ------------------------------------------%%
set(gcf,'Color','w')
drawnow

exportgraphics(gcf,'output.tif','Resolution',600);

disp('Spearman correlation panel exported successfully.');
