% Created by H.B. on 2016/06/21
% This function takes in a path and computes its derivative using the
% Savitzky-Golay interpolation technique
% x is the position
% t is the time

function [SG1,SG0] = getDerivativeOfPath(x,derivativeWindowLength,timeStepUnit)

% get first and second derivatives
N = 1;                 % Order of polynomial fit
F = derivativeWindowLength;                % Window length
dx = 1/timeStepUnit; % time-step for pixel size
[~,g] = sgolay(N,F);   % Calculate S-G coefficients    
HalfWin  = ((F+1)/2) -1;
queryPoints = (F+1)/2:length(x)-(F+1)/2; 
SG0 = zeros(length(queryPoints),1); 
SG1 = zeros(length(queryPoints),1);
for n = queryPoints
  % Zeroth derivative (smoothing only)
  SG0(n) = dot(g(:,1),x(n - HalfWin:n + HalfWin));
  % 1st differential
  SG1(n) = dot(g(:,2),x(n - HalfWin:n + HalfWin));
%   % 2nd differential
%   SG2(n) = 2*dot(g(:,3)',y(n - HalfWin:n + HalfWin))';
end
SG1 = SG1/dx;         % Turn differential into derivative
% SG2 = SG2/(dx*dx);    % and into 2nd derivative

end