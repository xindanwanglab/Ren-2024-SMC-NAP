% Created by H.B. on 2016/06/01
% This function takes in Hi-C data at 26 degrees and obtains the speed
% along each chromosomal arm

% ori-ter speed, going up the matrix: maximum intensity along the rows
close all; clear all; clc; 
%addpath('C:\Users\Hugo\Documents\GitHub\Harvard\Interactome\BrandaoMatLib\');
%addpath('C:\Users\Hugo\Documents\GitHub\Harvard\Interactome\');
addpath('/Users/zhongqingren/Library/CloudStorage/OneDrive-SharedLibraries-IndianaUniversity/Wang, Xindan - Wang Lab Common OneDrive/Lab Folders/WARREN (ZHONGQING)/Data/HiC Warren analysis/SMC_rates/matrices/');
addpath('/Users/zhongqingren/Library/CloudStorage/OneDrive-SharedLibraries-IndianaUniversity/Wang, Xindan - Wang Lab Common OneDrive/Lab Folders/WARREN (ZHONGQING)/Data/HiC Warren analysis/SMC_rates/BrandaoMatlib/');
strain = 'BWX5603_seq55_';
% plotting options
msz = 50;
plot10_or_5 = 10; 
maxColor = 0.01; 
printFig = 1;
tungColorMap = tungColorScheme();
% [255 255 255
%    162 192 222
%    140 137 187
%    140 87 167
%   140 45 143
%    120 20 120
%    90 15 90
%    60 10 60
%    30 5 30
%    0 0 0]; % get Tung's color map
% scales = [1.644853,1.96,2.575829,3.290527,3.89059,4.417173,4.891676 ]; % water-level scales 
% usePointOfMaxDistance = 1; 
% 
% % fitting options
% strelVal = 1;
scales = [0.26,0.51,0.76,1.01,1.26];%,1.51];%[0.51,0.76,1.01,1.26,1.51,1.96,2.575829,3.290527]%,3.89059,4.417173,4.891676]; % water-level scales 
usePointOfMaxDistance = 0; 
time = [15 25 35]';
time_extrapolate = [0,35];

% fitting options
useImclose = 1;
doSubtraction = 0;
% size for more stringent filtering; size of the diamond
strelVal = [0,1,5,10,12,12];

% genome values
genomeLength = 4033459/1e3;
parS = 0;%round((4033459*359/360-1.1204e4)/1e3);% -1 deg parS site location in kb

load(strcat(strain,'BlockZippingPath_10k.mat')); % 10 kb data 
numFiles=length(udata.FileNames); % for 10 kb
filesList = 1:numFiles; 

% 
% desiredFigurePath = strcat('C:\Users\Hugo\Dropbox (Personal)\Hugo+Xindan+David\',...
%     '201605 Hugo analysis for tethering paper\Figure 1\',...
%     'Analysis\Final_figures_',num2str(plot10_or_5),'k_curvatureCorrected\');
% desiredFigurePath = strcat('C:\Users\Hugo\Dropbox (Personal)\Hugo+Xindan+David\',...
%     '201605 Hugo analysis for tethering paper\Figure 1\',...
%     'Analysis\Plain_thresholding\');
desiredFigurePath = strcat('./',strain(1:end-1),'/');
numScales = length(scales); 
maxExtent = zeros((numFiles-1)*numScales,5); 
M = udata.HiC.matrix{1};
x = ((0:size(M,1))-size(M,1)/2)*udata.HiC.binsize/1e3; % in kb

% get maximum extent by subtracting adjacent images
% get maximum extent by water-level plot
% get maximum extent by subtracting from t=0
% ter-ori speed, going down the matrix: get maximum intensity along the
% columns


% loop through all files, all scales
count = 0;
for f = filesList 
    fc=1;
    C = udata.HiC.matrix{fc}; % reference "unzipped" data
    M = udata.HiC.matrix{f}; % time-point data
    
    if doSubtraction==1 
        R = (M-C); % subtract reference data
    else
        R = (M);
    end
    
    
    Rstd = zeros(size(R)); % standard error matrix (for diagonals)
    Rmean = zeros(size(R)); % expected contact matrix 
    
%     % compute standard error for every diagonal Hi-C map slice
%     L = size(R,1);
%     offset=0;
%     while offset<L-offset
%         offset= offset+1;
%         % get values of diagonal
%         vals=[diag(R,offset);diag(R,L-offset)];
%         % robust estimation of std and median
%         stdVal = 1.4826*mad(vals,1);
%         meanVal = nanmedian(vals);
%         Rstd=Rstd+diag(stdVal*ones(L-offset,1),offset);
%         Rmean=Rmean+diag(meanVal*ones(L-offset,1),offset);
%         
%         if L-offset>offset+1
%             Rstd=Rstd+diag(stdVal*ones(offset,1),L-offset);
%             Rmean=Rmean+diag(meanVal*ones(offset,1),L-offset);
%         end    
%     end

        vals = 1.4826*mad(R(:));
        meanVal = median(R(:));
        Rstd = vals*ones(size(R))/2;
        Rmean = meanVal*ones(size(R))/2;
        
    Rstd = Rstd+Rstd';
    Rmean = Rmean+Rmean';
    Rmean(Rmean==0) = max(Rmean(:));

    % compute water-level for each scale/water-level height
    for sc=scales 
        count = count+1;
        pval = 1-(normcdf(sc,0,1)-normcdf(-sc,0,1)); % p-value of given sigma 
        if doSubtraction==1
            R = fftshift(M-C); 
        else
            R = fftshift(M );
        end

        
        % obtain largest contiguous connected region
        R(R<=(Rmean+Rstd*sc)) =  0; % apply water-level cutoff
        R(R>0) =  1; 
        
        % mask out the center of R
        sR = size(R,1);
        mask = triu(ones(sR),msz)'+triu(ones(sR),msz);
        mask = mask.*(ones(sR)-triu(ones(sR),sR-msz)-triu(ones(sR),sR-msz)');
        R = R.*mask;
            
            
        if useImclose == 1
            se = strel('diamond',strelVal(scales==sc));
            R=imclose(R,se); 
        end        
        
                
        listOfPixels = growRegion(R); % connect/find largest region
        R(listOfPixels) = 2; 
        
        
        [I,J] = ind2sub(size(R),listOfPixels); 
        
        if usePointOfMaxDistance==1
%             % get point of furthest distance from parS site
%             parSbinX = size(M,1)/2+parS/udata.HiC.binsize*1e3; 
%             parSbinY = size(M,1)/2-parS/udata.HiC.binsize*1e3; 
%             [~,ind] = max((I-parSbinY).^2+(J-parSbinX).^2);
%             maxExtent(count,:) = [f,sc pval x(J(ind)),x(I(ind))];

            parSbin = size(M,1)/2+1+parS/udata.HiC.binsize*1e3; 
            lindists = abs(I-parSbin)+abs(J-parSbin);
            [val,~] = max(lindists); % furthest linear separation
            inds = find(lindists == val);
            maxExtent(count,:) = [f,sc pval median(x(J(inds))),median(x(I(inds)))];
        else
            % get point of furthest extent on x-axis and y-axis
            maxExtent(count,:) = [f,sc pval x(min(J)),x(max(I))];
        end
        
        % scale plots to desired color
        MM = M;
        MM(MM>maxColor) = maxColor; % set threshold color (top)
        MM(MM==0)= max(MM(:)); % set zeros (main diagonal) to max color
                
        % make plots from data
        h = figure(count);
        subplot(1,2,1);
        imagesc(x,x,(fftshift(MM))); axis square; set(gca,'ydir','normal');
        colormap(gcf,tungColorMap);
        xlabel('Genome position (kb)')
        ylabel('Genome position (kb)')
        
        % plot contiguous regions
        subplot(1,2,2);
        imagesc(x,x,(R));  set(gca,'ydir','normal'); axis square; 
        colormap(gca,'parula');
        xlabel('Geonomic position (kb)');
        ylabel('Geonomic position (kb)');
        title(strcat('p=',num2str(pval,1),' or \sigma=',num2str(sc,2)));

        % superimpose lines of maximum extent
        hold on; 
        plot(x,ones(size(x))*maxExtent(count,4),'k--','linewidth',1.5);
        plot(ones(size(x))*maxExtent(count,5),x,'k--','linewidth',1.5);
        plot(x,ones(size(x))*maxExtent(count,5),'k--','linewidth',1.5);
        plot(ones(size(x))*maxExtent(count,4),x,'k--','linewidth',1.5);
        set(gca,'fontsize',16);
        
        subplot(1,2,1);
        hold on;
        plot(x,ones(size(x))*maxExtent(count,4),'k--','linewidth',1.5);
        plot(ones(size(x))*maxExtent(count,5),x,'k--','linewidth',1.5);
        plot(x,ones(size(x))*maxExtent(count,5),'k--','linewidth',1.5);
        plot(ones(size(x))*maxExtent(count,4),x,'k--','linewidth',1.5);
        set(gca,'fontsize',16);
        if printFig==1 % save figure to pdf
            set(h,'PaperOrientation','landscape');
            set(h,'PaperUnits','normalized');
            set(h,'PaperPosition', [0 0 1 1]);
            print(h,strcat(desiredFigurePath,udata.FileNames{f},'_waterLevelPlot_sigma_',num2str(sc,3),'.pdf'),'-dpdf'); 
        end
    end
end

%%
% close all;
% plot distance vs time
% time = [30 40 50 ]';
for sc = scales
    pval = 1-(normcdf(sc,0,1)-normcdf(-sc,0,1)); % p-value for associated sigma

    h= figure(); % new figure
    ind = maxExtent(:,2)== sc; 
    
    terOri = min(maxExtent(ind,4:5)')'; 
    oriTer = max(maxExtent(ind,4:5)')';
    
    % plot maximum extent on each arm separately
    hh(1) = plot(time,terOri,'o','linewidth',2); hold on; 
    hh(2) =plot(time,oriTer,'^','linewidth',2); 
    
    % get speeds from linear fits (ter->ori direction)
    X = [ones(size(time)) time];
    t = time_extrapolate;
    [p,err] = lscov(X,terOri);
    hold on; plot(t,p(2)*(t)+(p(1)),'k','linewidth',1); 
    plot(t,(p(2)-err(2))*(t)+(p(1)+err(1)),'--k','linewidth',1);
    plot(t,(p(2)+err(2))*(t)+(p(1)-err(1)),'--k','linewidth',1);

    % get speeds from linear fits (ori->ter direction)
    X = [ones(size(time)) time];
    t = time_extrapolate;
    [pR,errR] = lscov(X,oriTer);
    hold on; plot(t,pR(2)*(t)+(pR(1)),'k','linewidth',1); 
    plot(t,(pR(2)-errR(2))*(t)+(pR(1)+errR(1)),'--k','linewidth',1);
    plot(t,(pR(2)+errR(2))*(t)+(pR(1)-errR(1)),'--k','linewidth',1);
    title(strcat('p=',num2str(pval,1),' or \sigma=',num2str(sc,2)));
    ylabel('Genomic position (kb)');
    xlabel('Time (min)');
    ylim([-2000,2000]);
    
    
    % parameters and errors
    param_list = {'Rate clockwise';'Rate counter-clockwise';...
                  'Rate error clockise';'Rate error counter-clockwise';...
                  'Lag-time clockwise';'Lag time counter-clockwise'; ...
                  'Lag-time error clockwise';'Lag-time error counter-clockwise'};
    lag = -(p(1)/p(2));
    lag_min = -(p(1)+err(1))/(p(2)-err(2));
    lag_max = -(p(1)-err(1))/(p(2)+err(2));
    lag_err = abs(lag_max-lag_min)/2;
    lagR = -(pR(1)/pR(2));
    lagR_min = -(pR(1)+errR(1))/(pR(2)-errR(2));
    lagR_max = -(pR(1)-errR(1))/(pR(2)+errR(2));
    lagR_err = abs(lagR_max-lagR_min)/2;
    param_list_values = [p(2),pR(2),err(2),errR(2),lag,lagR,lag_err,lagR_err]';
    %%%%%%%%%%%
    
%     legend(hh([2,1]),{strcat('Rate=',num2str(round(p(2))),'\pm',num2str(round(err(2))),'kb/min'),...
%     strcat('Rate=',num2str(round(pR(2))),'\pm',num2str(round(errR(2))),'kb/min')},...
%     'location','northwest');
    legend(hh([1,2]),{strcat('Rate=',num2str(round(p(2))),'\pm',...
        num2str(round(err(2))),'kb/min, ',...
        'Lag=',num2str(round(lag)),'\pm',...
        num2str(round(lag_err)),'min'),...
    strcat('Rate=',num2str(round(pR(2))),'\pm',num2str(round(errR(2))),'kb/min, ',...
        'Lag=',num2str(round(lagR)),'\pm',...
        num2str(round(lagR_err)),'min')},...
    'location','northwest');
    grid on;
    set(gca,'fontsize',16);
    if printFig==1
    print(h,strcat(desiredFigurePath,'DistanceVsTime_waterLevelPlot_sigma_',num2str(sc,3),'.pdf'),'-dpdf'); 
    end
    
    if printFig==1
        T = table(time,terOri,oriTer); 
        TT = table(param_list,param_list_values);
        writetable(T,strcat(desiredFigurePath,'posvstime_coordinates_',num2str(sc),'.txt'),'Delimiter','\t');
        writetable(T,strcat(desiredFigurePath,'parameter_fits_',num2str(sc),'.txt'),'Delimiter','\t');
    end
end
%%
if printFig==1
    clear h; 
    save(strcat(desiredFigurePath,'Analysis_parameters.mat')); 
end