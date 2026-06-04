%%====================================================================%%
% Flash Drought Onset Timing Coherence Panel
% ΔT (Onset Timing Difference) + OCR (Onset Coincidence Ratio)
%
% Author  : Vaibhav Kumar
% Project : Flash Drought India
% Created : 06 March 2026
%%====================================================================%%

close all; clear; clc;

fprintf('\nFlash Drought Onset Timing Coherence Panel\n');
fprintf('Run date: %s\n\n',datestr(now));

%%-----------------------Load data-------------------------%%
load('WATERGAP_vs_ERA5L_DeltaT.mat','DT');   DT_WG  = DT;
load('H08_vs_ERA5L_DeltaT.mat','DT');        DT_H08 = DT;
load('GLEAM_vs_ERA5L_DeltaT.mat','DT');      DT_GLE = DT;

load('WATERGAP_vs_ERA5L_OCR.mat','OCR');     OCR_WG  = OCR;
load('H08_vs_ERA5L_OCR.mat','OCR');          OCR_H08 = OCR;
load('GLEAM_vs_ERA5L_OCR.mat','OCR');        OCR_GLE = OCR;

%%---------------------Domain definition--------------------%%

lon0 = 65; lon1 = 96.5;
lat0 = 5;  lat1 = 36.5;

[nr,nc] = size(DT_WG);

lon = linspace(lon0,lon1,nc);
lat = linspace(lat0,lat1,nr);

%%--------------------------Grid style---------------------%%

degStep   = 5;
gridColor = [0.5 0.5 0.5];
gridAlpha = 0.20;

bgColor = [0.95 0.95 0.95];

applyGrid = @(ax) set(ax,...
'XGrid','on','YGrid','on',...
'GridColor',gridColor,'GridAlpha',gridAlpha,...
'XTick',lon0:degStep:lon1,...
'YTick',lat0:degStep:lat1,...
'TickDir','out',...
'Box','on',...
'FontSize',12,...
'FontWeight','bold');

%%-------------------------Color limits---------------------%%

CLIM_DT  = [0 2];
CLIM_OCR = [0 1];

%%----------------------Figure-----------------------------%%

figure('Units','centimeters','Position',[2 2 32 20])

t = tiledlayout(2,3,'TileSpacing','compact','Padding','loose');

colormap(parula(256))

DATA = {DT_WG,DT_H08,DT_GLE,...
        OCR_WG,OCR_H08,OCR_GLE};

titles = {'WATERGAP vs ERA5 Land',...
          'H08 vs ERA5 Land',...
          'GLEAM vs ERA5 Land'};

labels = {'(a1)','(a2)','(a3)','(b1)','(b2)','(b3)'};

%%----------------------Plot panels------------------------%%

for i = 1:6
    
    nexttile
    
    img = imagesc(lon,lat,DATA{i});
    set(img,'AlphaData',~isnan(DATA{i}));
    
    set(gca,'YDir','normal','Color',bgColor)
    
    axis image
    axis tight
    
    if i <= 3
        clim(CLIM_DT)
    else
        clim(CLIM_OCR)
    end
    
    set(gca,'LineWidth',1.2)
    
    ax = gca;
    applyGrid(ax)
    
    % Column titles only for first row
    if i <= 3
        title(titles{i},'FontWeight','bold')
    end
    
    % Panel label
    text(0.02,0.98,labels{i},'Units','normalized',...
        'FontSize',12,'FontWeight','bold','VerticalAlignment','top');
    
end

%%-----------------Row labels (ΔT and OCR)------------------%%

annotation('textbox',[0.02 0.68 0.04 0.05],...
'String','\DeltaT','EdgeColor','none',...
'FontSize',14,'FontWeight','bold','Rotation',90);

annotation('textbox',[0.02 0.30 0.04 0.05],...
'String','OCR','EdgeColor','none',...
'FontSize',14,'FontWeight','bold','Rotation',90);

%%-------------------ΔT colorbar (top row)------------------%%

drawnow

ax = findall(gcf,'Type','axes');
ax = flipud(ax);

ax3 = ax(3);
ax1 = ax(1);

pos3 = ax3.Position;
pos1 = ax1.Position;

bottom = pos3(2);
top    = pos1(2)+pos1(4);
height = top-bottom;

left   = pos3(1)+pos3(3)+0.02;
width  = 0.018;

cb1 = colorbar('Position',[left bottom width height]);
cb1.Ticks = 0:0.5:2;
cb1.FontSize = 12;
cb1.FontWeight = 'bold';
cb1.Label.String = '\DeltaT (pentads)';
cb1.Label.FontWeight = 'bold';

%%-----------------OCR colorbar (bottom row)----------------%%

ax6 = ax(6);
ax4 = ax(4);

pos6 = ax6.Position;
pos4 = ax4.Position;

bottom = pos6(2);
top    = pos4(2)+pos4(4);
height = top-bottom;

left   = pos6(1)+pos6(3)+0.02;

cb2 = colorbar('Position',[left bottom width height]);
cb2.Ticks = 0:0.2:1;
cb2.FontSize = 12;
cb2.FontWeight = 'bold';
cb2.Label.String = 'Onset coincidence ratio';
cb2.Label.FontWeight = 'bold';

%% =========================================================
%% Export figure
%% =========================================================

set(gcf,'Color','w')

exportgraphics(gcf,...
'RZSM_OnsetTimingCoherence_India_final.tif','Resolution',600)

disp('ΔT + OCR panel exported successfully.')
