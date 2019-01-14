function [fitresult, gof] = createFits(capacity, crate, resistance)
%CREATEFITS(CAPACITY,CRATE,RESISTANCE)
%  Create fits.
%
%  Data for 'untitled fit 1' fit:
%      X Input : capacity
%      Y Input : crate
%      Z Output: resistance
%  Output:
%      fitresult : a cell-array of fit objects representing the fits.
%      gof : structure array with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 09-Jan-2019 14:48:58

%% Initialization.

% Initialize arrays to store fits and goodness-of-fit.
fitresult = cell( 2, 1 );
gof = struct( 'sse', cell( 2, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

%% Fit: 'untitled fit 1'.
[xData, yData, zData] = prepareSurfaceData( capacity, crate, resistance );

% Set up fittype and options.
ft = fittype( 'k/(x*(a*x+b*y)^d)', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMinChange = 1e-10;
opts.Display = 'Off';
opts.MaxIter = 500;
opts.Robust = 'Bisquare';
opts.StartPoint = [0.626900804567717 0.936996636793824 0.514080993552528 0.124354043910244];

% Fit model to data.
[fitresult{1}, gof(1)] = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult{1}, [xData, yData], zData );
legend( h, 'untitled fit 1', 'resistance vs. capacity, crate', 'Location', 'NorthEast' );
% Label axes
xlabel capacity
ylabel crate
zlabel resistance
grid on
view( 1.7, 22.8 );

%% Fit: 'untitled fit 2'.
% Cannot generate code for fit 'untitled fit 2' because the data selection is incomplete.

