%% =========================================================
%% Multi-Product Consensus Probability for Flash Drought
%% FD_India Project
%% =========================================================

close all; clear; clc

fprintf('Starting consensus probability computation...\n')

%% =========================================================
%% Load Flash Drought Masks (1 = FD, 0 = No FD)
%% =========================================================

load('WATERGAP_FDmask_India.mat','FD_mask');  WG = FD_mask;
load('H08_FDmask_India.mat','FD_mask');       H08 = FD_mask;
load('ERA5L_FDmask_India.mat','FD_mask');     ERA = FD_mask;
load('GLEAM_FDmask_India.mat','FD_mask');     GLEAM = FD_mask;

clear FD_mask

%% =========================================================
%% Convert dimension: lat × lon × time → time × lat × lon
%% =========================================================

WG    = permute(WG,[3 1 2]);
H08   = permute(H08,[3 1 2]);
ERA   = permute(ERA,[3 1 2]);
GLEAM = permute(GLEAM,[3 1 2]);

%% =========================================================
%% Harmonize time dimension
%% =========================================================

T = min([size(WG,1), size(H08,1), size(ERA,1), size(GLEAM,1)]);

WG    = WG(1:T,:,:);
H08   = H08(1:T,:,:);
ERA   = ERA(1:T,:,:);
GLEAM = GLEAM(1:T,:,:);

[nt,ny,nx] = size(ERA);

fprintf('Final harmonized dataset: %d time × %d lat × %d lon\n',nt,ny,nx);

%% =========================================================
%% Correct latitude orientation (north-up)
%% =========================================================

WG    = WG(:,end:-1:1,:);
H08   = H08(:,end:-1:1,:);
ERA   = ERA(:,end:-1:1,:);
GLEAM = GLEAM(:,end:-1:1,:);

%% =========================================================
%% Ensure binary mask (safety step)
%% =========================================================

WG    = double(WG>0);
H08   = double(H08>0);
ERA   = double(ERA>0);
GLEAM = double(GLEAM>0);

%% =========================================================
%% Stack products
%% dimension → time × lat × lon × product
%% =========================================================

FD_stack = cat(4,WG,H08,ERA,GLEAM);

%% =========================================================
%% Compute consensus probability
%% Equation: P_cons(i,t) = (1/4) Σ I_k
%% =========================================================

P_cons = mean(FD_stack,4,'omitnan');

%% =========================================================
%% High-confidence flash drought detection
%% ≥ 3 products agree
%% =========================================================

FD_high_conf = P_cons >= 0.75;

%% =========================================================
%% Time-averaged consensus probability (final spatial map)
%% =========================================================

P_cons_mean = squeeze(mean(P_cons,1,'omitnan'));

%% =========================================================
%% Diagnostics
%% =========================================================

fprintf('Consensus probability range: %.2f – %.2f\n', ...
        min(P_cons(:)), max(P_cons(:)));

fprintf('Mean consensus range: %.2f – %.2f\n', ...
        min(P_cons_mean(:)), max(P_cons_mean(:)));

%% =========================================================
%% Save outputs
%% =========================================================

save('FD_ConsensusProbability_India.mat',...
     'P_cons','FD_high_conf','P_cons_mean','-v7.3')

fprintf('Consensus probability successfully computed and saved.\n')

%% =========================================================
%% Quick visualization (optional check)
%% =========================================================

figure
imagesc(P_cons_mean)
axis image
axis off
colormap(parula(256))
colorbar
title('Mean Consensus Probability')