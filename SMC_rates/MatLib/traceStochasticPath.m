% Created by H.B. on 01/09/2015
% This function takes in a HiC matrix of values and starting point. It then
% traces a path of maximum likelihood of set length (using nearest
% neighbours only)
% M is the matrix, x is the initial point

function [meanPath_X,meanPath_Y, paths] = traceStochasticPath(M,numPaths,rateOfScramble,makePlot)
    maxL = size(M,1); 
    if ~exist('L','var'); 
       L = maxL;  
    end
    if ~exist('numPaths','var'); 
       numPaths = 1;  
    end
    if ~exist('rateOfScramble','var'); 
       rateOfScramble = 0.1;  
    end
    if ~exist('makePlot','var'); 
       makePlot = 0;  
    end
    if makePlot~=0 
        h = figure(); imagesc(M); set(gca,'ydir','normal'); axis square; 
    end
    paths = cell(numPaths,1);
    % stop when any of the coordinates L-n, n are reached
    for p=1:numPaths
    clear path; 
    pos = [1,1]; 
    idx = 0; 
    while (L-pos(2))>=pos(1) & pos(2)<L & pos(1)<L 
        idx = idx+1;
        % get nearest neighbour
        newPos = [pos(1),pos(2)+1; % right
                  pos(1)+1,pos(2); % up
                  pos(1)+1,pos(2)+1]; % up-right
        vals = M(sub2ind(size(M),newPos(:,1),newPos(:,2))); 
        vals(vals==0) = eps;
        cdfs = cumsum(vals./sum(vals)); 
        
        % get new position 
        if rand()<rateOfScramble
            [~,dxn] = max(cdfs>=rand()); % add some noise to path
        else
            [~,dxn] = max(vals./sum(vals));
        end         
        path(idx) = {pos};
        pos = newPos(dxn,:); 
    end
    paths{p}=cell2mat(path');
    if makePlot~=0
        figure(h); 
        hold on; 
        plot(paths{p}(:,2),paths{p}(:,1),'k'); 
    end
    end  
    % get mean paths
    paths2 = [];
    for p=1:numPaths
        paths2 = [paths2; [ paths{p}(:,1),paths{p}(:,2)]];%#ok<*AGROW> %-paths{p}(:,1)]];
    end
    minXrange  = min(paths2(:,1));
    maxXrange = max(paths2(:,1));
    minYrange  = min(paths2(:,2));
    maxYrange = max(paths2(:,2));   
 
    meanPath_X = zeros(maxXrange-minXrange+1,1);
    stdPath_X = meanPath_X; 
    for p=1:(maxXrange-minXrange+1)
        inds = paths2(:,1)==p;
        meanPath_X(p) = mean(paths2(inds,2));
        stdPath_X(p) = std(paths2(inds,2));%/sqrt(sum(inds));
    end
    meanPath_Y = zeros(maxYrange-minYrange+1,1);
    stdPath_Y = meanPath_Y; 
    for p=1:(maxYrange-minYrange+1)
        inds = paths2(:,2)==p;
        meanPath_Y(p) = mean(paths2(inds,2));
        stdPath_Y(p) = std(paths2(inds,2));%/sqrt(sum(inds));
    end 
end