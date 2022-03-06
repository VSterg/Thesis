function plotmts(xM)
% Plots multivariate time series in one figure using subplots
height = size(xM,2); % number of variables stored in height
for k = 1:height % make a plot for each variable
    subplot(height,1,k)
    plot(xM(:,k))
end
