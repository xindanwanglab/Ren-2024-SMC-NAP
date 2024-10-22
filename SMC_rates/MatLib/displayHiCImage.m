% Created by H.B. on 2016/06/14
% Takes in a single image or series of images as a cell array
% default color map is Tung's color map
function displayHiCImage(img,xy_label,img_title,xy_ticks,ax_handles)
fig = gcf(); 

if ~exist('xy_ticks','var')
    xy_ticks = 1:size(img,1); 
end

stopScale = 3; 
tungColorMap = tungColorScheme();

S = img; 
stdVal = nanstd(S(:));
maxVal = nanmedian(S(:)) + stdVal*stopScale;

S(S>maxVal) = maxVal; 
if exist('ax_handles','var')
    fig.CurrentAxes = ax_handles;
end
imagesc(xy_ticks,xy_ticks,S); axis square;set(gca,'ydir','normal'); 
xlabel(xy_label);
ylabel(xy_label);
title(img_title,'interpreter','none');
set(gca,'fontsize',16); 
colormap(fig,tungColorMap);
end