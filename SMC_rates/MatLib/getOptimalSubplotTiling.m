% Created by H.B. on 15/07/2015
% This function takes as arguments the total number of plots to be made and
% based on the screen resolution calculates the optimal tiling of subplots
% for ease of visualization via minimization of unused area

function [tiling,M] = getOptimalSubplotTiling(numPlots,fig)
    % get figure dimensions
%     X = get(0,'screensize');
    X = get(fig,'position');
    dimX = X(3); 
    dimY = X(4); 
    fitness = inf*ones(numPlots,numPlots); 
    for sx = 1:numPlots
       for sy = 1:numPlots
           sqd = min(dimY/sy,dimX/sx);
           if (sx*sy >= numPlots)
               % minimize left-over area, maximize used area
               fitness(sy,sx) = (dimX*dimY-sx*sy*sqd^2) + (sx*sy-numPlots)*sqd^2;     
           end
       end
    end
    [M,I] = min(fitness(:));
    tiling = zeros(1,1); 
    [tiling(1), tiling(2)] = ind2sub(size(fitness),I);
end


