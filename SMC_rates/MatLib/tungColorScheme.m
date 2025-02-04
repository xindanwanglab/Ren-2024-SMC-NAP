function tungColorMap = tungColorScheme();
map =[255, 255, 255;...
    162, 192, 222;...
    140, 137, 187;...
    140, 87, 167;...
    140, 45, 143;...
    120, 20, 120;...
    90, 15, 90;...
    60, 10, 60;...
    30, 5, 30;...
    0, 0, 0]./255;

x = (0:7.11:64)';
xq = (0:1+1/64:64)';

tungColorMap = zeros(64,3);

for dim=1:3;
tungColorMap(:,dim) = interp1(x,map(:,dim),xq);
end;
return; 