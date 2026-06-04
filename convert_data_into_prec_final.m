%%------------------------------------------------------------%%
% FD_India: Pixel-wise Pentad-of-Year Percentile Transformation
%%------------------------------------------------------------%%

clear; close all; clc;

%%---------------- USER SETTINGS ----------------%%
inFile  = "file01.mat";  
varName = "RZSM";                                   
minValidYears = 5;   % Recommended: 20 (out of 40)

%%---------------- LOAD -------------------------%%
S = load(inFile);
X = double(S.(varName));

X(~isfinite(X)) = NaN;
X(X < 0) = NaN;

[nRows, nCols, nTime] = size(X);
nPentadsPerYear = 73;

assert(nTime == 73*40, 'Time dimension must equal 2920 (73×40).');

fprintf("Data size: %d x %d x %d\n", nRows, nCols, nTime);

%%---------------- PENTAD INDEX ----------------%%
pentad_of_year = mod((1:nTime)-1, nPentadsPerYear) + 1;

%%---------------- OUTPUT INIT -----------------%%
RZSM_pct = NaN(nRows, nCols, nTime);

%%---------------- CORE COMPUTATION ------------%%
for p = 1:nPentadsPerYear
    
    tIdx = find(pentad_of_year == p);      
    V = X(:,:,tIdx);                       
    
    V2 = reshape(V, [], length(tIdx));     % pixels × years
    
    validMask = ~isnan(V2);
    nValid = sum(validMask, 2);
    
    [Vsort, sortIdx] = sort(V2, 2, 'ascend', 'MissingPlacement','last');
    
    P = NaN(size(V2));
    ok = nValid >= minValidYears;
    
    for i = find(ok)'
        
        ni = nValid(i);
        if ni == 0
            continue
        end
        
        % Bias-corrected empirical percentile
        ri = 1:ni;
        pi = 100 * ((ri - 0.5) / ni);   % <-- FIXED HERE
        
        P(i, sortIdx(i,1:ni)) = pi;
    end
    
    RZSM_pct(:,:,tIdx) = reshape(P, nRows, nCols, length(tIdx));
    
    fprintf("Pentad %d / 73 completed\n", p);
end

%%---------------- SAVE ------------------------%%
save("ERA5L_RZSM_1981-2020_India_Pentad_Percentile.mat","RZSM_pct", "-v7.3");
disp("Pentad-of-year percentile transformation completed.");

%%---------------- Percentile Distribution Check ----------------%%
% Flatten all valid percentile values
P_all = RZSM_pct(:);
P_all = P_all(~isnan(P_all));

fprintf('Min percentile: %.3f\n', min(P_all));
fprintf('Max percentile: %.3f\n', max(P_all));
fprintf('Mean percentile: %.3f\n', mean(P_all));
fprintf('Std percentile: %.3f\n', std(P_all));

%%---------------- Histogram ----------------%%
figure;
histogram(P_all, 50, 'Normalization','pdf');
xlabel('Percentile');
ylabel('Probability Density');
title('Distribution of RZSM Percentiles (1981–2020)');
grid on;
