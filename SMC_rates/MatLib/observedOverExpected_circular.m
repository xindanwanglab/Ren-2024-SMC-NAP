% Created by H.B. on 2016/06/14
% This function calculates the observed over expected Hi-C map 
% R is the input image

function [R,Rmean,Rstd] = observedOverExpected_circular(R)
Rstd = zeros(size(R));
Rmean = zeros(size(R));
L = size(R,1);
offset=0;
while offset<L-offset
    offset= offset+1;
    vals=[diag(R,offset);diag(R,L-offset)];
    stdVal = nanstd(vals);
    meanVal = nanmean(vals);
    Rstd=Rstd+diag(stdVal*ones(L-offset,1),offset);
    Rmean=Rmean+diag(meanVal*ones(L-offset,1),offset);
    if L-offset>offset+1
    Rstd=Rstd+diag(stdVal*ones(offset,1),L-offset);
    Rmean=Rmean+diag(meanVal*ones(offset,1),L-offset);
    end    
end
Rstd = Rstd+Rstd';
Rmean = Rmean+Rmean';
Rmean(Rmean==0) = max(Rmean(:));
R = R./Rmean;     
end