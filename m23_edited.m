    % Flash drought detection using MATLAB
    % Vaibhav Kumar version 2    2024/3/23
    %         version 3     2024/4/03
    %         version 3.1   2024/4/19
    %         version 3.2   2024/4/22
    % source from https://github.com/Hydroclimate2023/global-flash-drought/tree/master/code
    % n, number of pentads 
    % a(n), pentad soil moisture percentile [0-100] 
    % b(1), number of flash droughts 
    % b(2), mean duration of flash droughts
    % b(3), mean severity of flash droughts   
    %       the mean accumulated soil moisture percentile deficits from the threshold, 40%
    % b(4), mean speed of flash droughts (percentile/pentads)
     
    clear all; close all; clc;
   
    load ("file.mat");
    
    bb1(1:size(b5,1),1:size(b5,2))=0;
    bb2(1:size(b5,1),1:size(b5,2))=0;
    bb3(1:size(b5,1),1:size(b5,2))=0;
    bb4(1:size(b5,1),1:size(b5,2))=0;
    tot_drt(1:size(b5,1),1:size(b5,2),1:size(b5,3))=0;

        for l=1:size(b5,1)
        for k=1:size(b5,2)
    
    a1(:)=b5(l,k,1:size(b5,3));
    if isnan(a1(1))  
        bb1(l,k)=NaN;
        bb2(l,k)=NaN;
        bb3(l,k)=NaN;
        bb4(l,k)=NaN;
        continue; 
    end
 
% % Data transform to percentile 
% %     for i=1:550
% %     value=a1(i);
% %     perc = prctile(a1,1:100);
% %     [c index] = min(abs(perc'-value));
% %     x(i) = index+1;
% %     end
    a=a1;   
    thresh  =40;     % drought start threshold [percentile]
    thresh1 =20;     % drought end threshold   [percentile]
    speed   =5;      % drought speed threshold [percentile]
    td = 4;          % minimum drought duration threshold [pentads]
    td1 = 12;        % maximum drought duration threshold [pentads]
    n = length(a);
    cnt = 0; flag = 0; si = 0;
    odur = zeros(1, n); osev = zeros(1, n); ospd = zeros(1, n);
    tmin = 99;
    flag1 = 0; flag2 = 0;
    drt = zeros(1, n);
    
    for i = 1:n
        if a(i) < 0
            error('Error: Input data cannot be negative.');
        end
    end
    
    for j = 2:n
        if a(j) < thresh
            if flag == 0 % drought start
                flag1 = 1;
                if j > 1
                    if a(j - 1) < thresh
                        flag1 = 0;
                    end
                end
                if flag1 == 1
                    flag = 1;
                    flag2 = 1; % onset period
                    cnt = cnt + 1;
                    si = j; % si is drought start time
                    osev(cnt) = thresh - a(j);
                    odur(cnt) = 1;
                    ospd(cnt) = thresh - a(j);
                    drt(j) = 1;
                    if a(j) < tmin
                        tmin = a(j);
                    end
                end
            else
                if flag2 == 1 % onset period
                    if thresh - a(j) >= (j - si + 1) * speed
                        if tmin <= thresh1
                            if a(j) >= a(j - 1) % stop onset period
                                flag2 = 0;
                                if a(j) <= thresh1 % enter recovery period
                                    osev(cnt) = osev(cnt) + (thresh - a(j));
                                    odur(cnt) = odur(cnt) + 1;
                                    drt(j) = 2;
                                    if a(j) < tmin
                                        tmin = a(j);
                                    end
                                else % no recovery period
                                    flag = 0;
                                    tmin = 99;
                                    if (odur(cnt) < td)|(odur(cnt) > td1)
                                        osev(cnt) = 0;
                                        odur(cnt) = 0;
                                        ospd(cnt) = 0;
                                        drt(si:j - 1) = 0;
                                        cnt = cnt - 1;
                                    end
                                end
                            else % continue onset period
                                osev(cnt) = osev(cnt) + (thresh - a(j));
                                odur(cnt) = odur(cnt) + 1;
                                ospd(cnt) = (thresh - a(j)) / (j - si + 1);
                                drt(j) = 1;
                                if a(j) < tmin
                                    tmin = a(j);
                                end
                            end
                        else % before entering thresh1, continue onset period
                            osev(cnt) = osev(cnt) + (thresh - a(j));
                            odur(cnt) = odur(cnt) + 1;
                            ospd(cnt) = (thresh - a(j)) / (j - si + 1);
                            drt(j) = 1;
                            if a(j) < tmin
                                tmin = a(j);
                            end
                        end
                    else % speed does not meet
                        flag2 = 0; % so enter recovery period
                        if tmin <= thresh1
                            if a(j) <= thresh1
                                osev(cnt) = osev(cnt) + (thresh - a(j));
                                odur(cnt) = odur(cnt) + 1;
                                drt(j) = 2;
                            else % drought dismiss
                                flag = 0;
                                tmin = 99;
                                if  (odur(cnt) < td)|(odur(cnt) > td1)
                                    osev(cnt) = 0;
                                    odur(cnt) = 0;
                                    ospd(cnt) = 0;
                                    drt(si:j - 1) = 0;
                                    cnt = cnt - 1;
                                end
                            end
                        else % tmin > thresh1, does not meet drought criterion
                            drt(si:j - 1) = 0;
                            osev(cnt) = 0;
                            odur(cnt) = 0;
                            ospd(cnt) = 0;
                            cnt = cnt - 1;
                            flag = 0;
                            tmin = 99;
                        end
                    end
                else % flag2 = 0, recovery period
                    if a(j) <= thresh1
                        osev(cnt) = osev(cnt) + (thresh - a(j));
                        odur(cnt) = odur(cnt) + 1;
                        drt(j) = 2;
                    else % drought dismiss
                        if  (odur(cnt) < td)|(odur(cnt) > td1)
                            osev(cnt) = 0;
                            odur(cnt) = 0;
                            ospd(cnt) = 0;
                            drt(si:j - 1) = 0;
                            cnt = cnt - 1;
                        end
                        flag = 0;
                        tmin = 99;
                    end
                end
            end
            if j == n && flag == 1
                if odur(cnt) < td || tmin > thresh1
                    osev(cnt) = 0;
                    odur(cnt) = 0;
                    ospd(cnt) = 0;
                    drt(si:j) = 0;
                    cnt = cnt - 1;
                end
            end
        else % a(j) >= thresh
            if flag == 1 % drought dismiss
                flag = 0;
                if odur(cnt) < td || tmin > thresh1
                    drt(si:j) = 0;
                    osev(cnt) = 0;
                    odur(cnt) = 0;
                    cnt = cnt - 1;
                end
                tmin = 99;
            end
        end
    end
    % odur, total duration of each flash drought
    % osev, total severity of each flash drought
    % ospd, onset speed of each flash drought
    if cnt > 0
        b(1) = cnt;
        b(2) = sum(odur(1:cnt)) / cnt;
        b(3) = sum(osev(1:cnt)) / cnt;
        b(4) = sum(ospd(1:cnt)) / cnt;
    else
        b = zeros(1, 4);
    end
    % output
    bb1(l,k)=b(1);
    bb2(l,k)=b(2);
    bb3(l,k)=b(3);
    bb4(l,k)=b(4);
    
    % drt, 1-drought onset, 2-drought recovery, 0-no drought
     tot_drt(l,k,1:size(b5,3))=drt;
end
    end
  
 %%-----------------Define target directory------------------%%
outDir = 'D:\FD_India\ERA5L_RZSM\ERA5L_FD_Characteristics';

%%----------------Save each variable separately in .mat format------------%%
save(fullfile(outDir,'ERA5L_FD_frequency_India.mat'),'bb1');
save(fullfile(outDir,'ERA5L_FD_mean-duration_India.mat'),'bb2');
save(fullfile(outDir,'ERA5L_FD_mean-severity_India.mat'),'bb3');
save(fullfile(outDir,'ERA5L_FD_mean-speed_India.mat'),'bb4');
save(fullfile(outDir,'ERA5L_FD_tot_drt_(drought-phase-flag)_India.mat'),'tot_drt');

%%----------------Optional confirmation----------------------------------%%
disp(['All output files saved successfully in: ', outDir]);

 figure; 
 h=imagesc(bb1);set(h,'alphadata',~isnan(bb1));axis off;axis equal;colorbar;colormap(flipud(hot(150)));%clim([0 40]);
 figure;
 h=imagesc(bb2);set(h,'alphadata',~isnan(bb2));axis off;axis equal;colorbar;colormap(flipud(hot(150)));%clim([0 30]);
 figure;
 h=imagesc(bb3);set(h,'alphadata',~isnan(bb3));axis off;axis equal;colorbar;colormap(flipud(hot(150)));%clim([0 600]);
 figure;
 h=imagesc(bb4);set(h,'alphadata',~isnan(bb4));axis off;axis equal;colorbar;colormap(flipud(hot(150)));%clim([0 40]);  

% figure;% 3 D view
% [xx, yy, zz] = meshgrid(1:135,1:129,1:550);
% xslice=200;yslice=200;
% zslice=[100, 200,300,400,500];
% h=slice(xx,yy,zz,tot_drt,xslice,yslice,zslice);
% set (h,'EdgeColor','none');axis equal; 
% set(gca,'YDir','reverse');

%%--------show time series for one locaton (no. of counts of FD events) (i.e. lon:80,lat:26)------------%%
%  figure;
% %%----------drt, 1-drought onset, 2-drought recovery, 0-nodrought-------------%%
%  dd(1:1050)=tot_drt(80,26,1:1050);
%  plot(1:1050, dd);
% 
% %%-------time series for input variable (i.e. focus on one location lon:80,lat:26)-------%%
%  figure;  
%  dd(1:1050)=b5(80,26,1:1050);
%  plot(1:1050, dd);
% 
% %%----------percentile distribution of input varibale (i.e. soil-moisture, evaopration, rainfall) on one time-----%% 
%  figure; 
%  h=imagesc(b5(:,:,33));set(h,'alphadata',~isnan(bb4));axis off;axis equal;colorbar;colormap(flipud(hot(150)));%clim([0 100]);
