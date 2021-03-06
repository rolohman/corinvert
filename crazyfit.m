function [fitresult, gof] = crazyfit(intd, intbp, hp, weights)
%CREATEFIT(INTD,INTBP,HP,WEIGHTS)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : intd
%      Y Input : intbp
%      Z Output: hp
%      Weights : weights
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 09-May-2019 14:46:43


%% Fit: 'untitled fit 1'.
[xData, yData, zData, weights_1] = prepareSurfaceData( intd, intbp, hp, weights );

% Set up fittype and options.
ft = fittype( 'a+(b+c*sin(x/360*2*pi-d))*y', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0 0 0 0];
opts.Weights = weights_1;

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );
return
% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, [xData, yData], zData );
legend( h, 'untitled fit 1', 'hp vs. intd, intbp with weights', 'Location', 'NorthEast' );
% Label axes
xlabel intd
ylabel intbp
zlabel hp
grid on
view( 21.6, 51.8 );


