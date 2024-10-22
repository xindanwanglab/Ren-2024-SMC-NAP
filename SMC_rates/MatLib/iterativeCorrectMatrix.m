function W = iterativeCorrectMatrix(M,plotTrue)
    if(~exist('plotTrue','var'))
        plotTrue=1; 
    end
        
    W = M;
    B = ones(1,length(W));
    for i=1:10
        S = nansum(W);
        dB = S./nanmean(S);
        B = B.*dB;
        W = (W./(dB'*dB));
%         W(~isfinite(W)) = 1;%max(M(:));
        %figure(i); plot(B); 
    end
    
    if plotTrue
        figure(); h=subplot(2,1,1); imagesc(M); set(gca,'ydir','normal');axis square; 
        hh=subplot(2,1,2); imagesc((W)); set(gca,'ydir','normal');axis square; 
        linkaxes([h hh])
%         figure(); plot(B);
    end
end