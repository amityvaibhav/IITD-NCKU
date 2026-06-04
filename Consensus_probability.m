%%====================================================================%%
% FD_India: Multi-Product Consensus Probability
%
% Author  : Prajwal Giri
% Project : FD_India
%%====================================================================%%

close all; clear; clc

fprintf('\nFD_India – Consensus Probability Analysis\n');
fprintf('Run date: %s\n\n',datestr(now));

%% =========================================================
%% Load consensus dataset
%% =========================================================

load('FD_ConsensusProbability_India.mat','P_cons','FD_high_conf','P_cons_mean');

%% =========================================================
%% Compute additional probability diagnostics
%% =========================================================

P_high_mean  = squeeze(mean(FD_high_conf,1,'omitnan'));   % ≥3 products
P_cons_event = squeeze(mean(P_cons >= 0.5,1,'omitnan'));  % ≥2 products

lat = size(P_cons_mean,1);
lon = size(P_cons_mean,2);

%% =========================================================
%% Mask extremely small probabilities
%% =========================================================

maskThr = 0.005;

P_cons_mean(P_cons_mean < maskThr)   = NaN;
P_high_mean(P_high_mean < maskThr)   = NaN;
P_cons_event(P_cons_event < maskThr) = NaN;

%% =========================================================
%% Save probability outputs
%% =========================================================

save('FD_MeanConsensusProbability_India.mat','P_cons_mean','-v7.3');
save('FD_HighConfidenceConsensusProbability_India.mat','P_high_mean','-v7.3');
save('FD_ModerateConsensusProbability_India.mat','P_cons_event','-v7.3');

fprintf('Probability maps saved successfully\n');

%% =========================================================
%% Determine color scaling
%% =========================================================

allVals = [P_cons_mean(:);P_high_mean(:);P_cons_event(:)];
allVals = allVals(~isnan(allVals));

CLIM = [0 prctile(allVals,98)];

%% =========================================================
%% Map coordinates
%% =========================================================

lon0 = 65; lon1 = 96.5;
lat0 = 5;  lat1 = 36.5;

lonVec = linspace(lon0,lon1,lon);
latVec = linspace(lat0,lat1,lat);

degStep   = 5;
gridColor = [0.6 0.6 0.6];
gridAlpha = 0.15;

%% =========================================================
%% Figure layout
%% =========================================================

figure('Units','centimeters','Position',[2 2 30 10])
set(gcf,'Color','w')   % fix figure background

t = tiledlayout(1,3,'TileSpacing','compact','Padding','loose');

colormap(turbo(256))

applyGrid = @(ax) set(ax,...
    'XGrid','on','YGrid','on',...
    'GridColor',gridColor,...
    'GridAlpha',gridAlpha,...
    'XTick',lon0:degStep:lon1,...
    'YTick',lat0:degStep:lat1,...
    'TickDir','out',...
    'Box','on',...
    'FontSize',11,...
    'FontWeight','bold');

titles = {'Mean Consensus Probability',...
          'High-Confidence Consensus (≥3 products)',...
          'Moderate Consensus (≥2 products)'};

labels = {'(a)','(b)','(c)'};

DATA = {P_cons_mean,P_high_mean,P_cons_event};

%% =========================================================
%% Plot panels
%% =========================================================

for i = 1:3

    nexttile

    img = imagesc(lonVec,latVec,DATA{i});
    set(img,'AlphaData',~isnan(DATA{i}));

    set(gca,'YDir','normal')
    set(gca,'Color','w')   % fix axes background

    axis image
    axis tight
    clim(CLIM)

    set(gca,'LineWidth',1.2)

    ax = gca;
    applyGrid(ax);

    title(titles{i},'FontWeight','bold')

    text(0.02,0.98,labels{i},...
        'Units','normalized',...
        'FontSize',12,...
        'FontWeight','bold',...
        'VerticalAlignment','top');

end

%% =========================================================
%% Shared vertical colorbar
%% =========================================================

drawnow

ax = findall(gcf,'Type','axes');
ax = flipud(ax);

posR = ax(3).Position;
posL = ax(1).Position;

bottom = posL(2);
top    = posL(2)+posL(4);
height = top-bottom;

left   = posR(1)+posR(3)+0.02;
width  = 0.018;

cb = colorbar('Position',[left bottom width height]);

cb.Ticks = linspace(CLIM(1),CLIM(2),6);
cb.FontSize = 13;
cb.FontWeight = 'bold';

cb.Label.String = 'Consensus Probability';
cb.Label.FontWeight = 'bold';
cb.Label.FontSize = 14;

%% =========================================================
%% Export figure
%% =========================================================

exportgraphics(gcf,...
'RZSM_ConsensusProbability_India_final.tif',...
'Resolution',600)

fprintf('Consensus probability panel exported successfully\n');