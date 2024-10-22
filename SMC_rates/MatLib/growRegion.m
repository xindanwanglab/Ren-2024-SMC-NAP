% Created by H.B. on 2016/06/01
% This function finds all connected points starting from an initial value

function listOfPixels = growRegion(img)%,startPoint)
%     indToLookFor = sub2ind(size(img),startPoint(2),startPoint(1));
    CC = bwconncomp(img);
%     listOfPixels = []; 
%     for n=1:CC.NumObjects
%         if sum(CC.PixelIdxList{n}==indToLookFor)>0
%             listOfPixels =  CC.PixelIdxList{n};
%             return; 
%         end
%     end
    
    % sort by size
    maxSize=0; 
    maxSizeIdx = 1; 
    for n=1:CC.NumObjects
        if length(CC.PixelIdxList{n})> maxSize
            maxSize = length(CC.PixelIdxList{n});
            maxSizeIdx = n; 
        end
    end
    listOfPixels = CC.PixelIdxList{maxSizeIdx}; 
end




