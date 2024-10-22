% Created by H.B. on 01/09/2015
% This function takes in a HiC matrix of values and starting point. It then
% traces a path of maximum likelihood of set length (using nearest
% neighbours only)
% M is the matrix, x is the initial point

function maxPathCoords = traceMinimumPathNN(M,x,L)
    maxL = size(M,1); 
    if ~exist('L','var'); 
       L = maxL;  
    end
    
    maxPathCoords = zeros(L,3); % [x y value]
    maxPathCoords(1,:) = [x,M(sub2ind(size(M),x(1),x(2)))];  
    pos = x; 
   
    % for the length of the path
    % always move left, up, or up-left to avoid looping back
    for p=2:L       
        % get nearest neighbour
        newPos = [pos(1),pos(2)-1; % left 
                  pos(1)+1,pos(2); % up
                  pos(1)+1,pos(2)-1]; % up-left 
%         newPos = [pos(1),pos(2)+1; % right
%                 pos(1)-1,pos(2); % down
%                 pos(1)-1,pos(2)+1]; % down-right
        % apply periodic boundary to these positions
        newPos = mod(newPos,maxL); 
        newPos(newPos == 0) = maxL; 
        
        % get nearest neighbour values and find maximum
        [val,ind] = min(M(sub2ind(size(M),newPos(:,1),newPos(:,2))));
        pos = newPos(ind,:); 
             
        % record position and value of nearest neighbour maximum
        maxPathCoords(p,:) = [pos, val];
    end
end