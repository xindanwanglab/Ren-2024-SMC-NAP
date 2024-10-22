x = ZipupratessummaryZRS5.timePoints;
Y = ZipupratessummaryZRS5{1:3,2:9};

figure; % Open a new figure window
hold on; % Hold on to the current plot

colors = ['b', 'r', 'g', 'k', 'c', 'm', 'y', 'b']; % Define a color for each set
markers = ['o', 'x', '+', '*', 's', 'd', '>', '^']; % Define a marker for each set
lines = '-'; % Define line styles
%names = ['BWX4070_R', 'BWX5062_R', 'BWX5063_R', 'BWX5663_R', 'BWX4070_L', 'BWX5062_L', 'BWX5063_L', 'BWX5663_L'];

for i = 1:size(Y,2)
    y = Y(:, i); % Extract the ith set of y data
    scatter(x, y, 100,strcat(colors(i), markers(i))); % Plot x vs. y with specific marker and color
    
    % Calculate and plot trend line
    p = polyfit(x, y, 1); % Fit a 1st-degree polynomial
    yfit = polyval(p, x); % Evaluate the polynomial at the original x values
    plot(x, yfit, strcat(colors(i), lines(1))); % Plot the trend line with specific color and line style
end

% Enhance plot
xlabel('Time(min)');
ylabel('Distance(kb)');
title('SMC translocation rates');
legendStrings = names;
legend(legendStrings{:}); % Generate dynamic legend based on the number of data sets
hold off; % Release the plot
