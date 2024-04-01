function [fit_result, y_fit, gof] = ...
    FUNC_noncentered_Lorentzian_fit_v1(x, y, w_guess, y_max_guess, c_guess)
%% DESCRIPTION: 
% Fits a Lorentzian known to be centered at 0. Fits for an overall scale factor y_max and a FWHM w
%
%% HISTORY:
% -v1 created on 8/31/19 by Dan Palken. Based directly upon FUNC_centered_Lorentzian_fit_v1
% -v2 created on 9/3/19 by DP. Fixes height at 1. 
%
%% VARIABLE INFO:
% INPUT
% - x: independant variable data. Uncertainties assumed appropriately low
% - y: dependant variable data
% - w_guess: inintial guess for the fit parameter w, corresponding to the full width at half max 
% (FWHM) of the Lorentzian
% - c_guess: initial guess for the parameter c, corresponding to the center frequency
% - weights (corresponding to inverse varaiances of the power gains) to be used in the fit
% OUTPUT
% - fit_result: a fit object representing the fit
% - y_fit: fitted values at all the input locations
% - gof: structure with goodness-of fit info
% - chi2_per_DOF: chi-squared per degree of freedom for the fit
%

%% ============================================================================================== %%
%% PREPARE FIT
% assumes that the center of the Lorentzian is at 0
[x_data, y_data] = prepareCurveData(x, y);
ft = fittype('w^2 / ((2 * (x-c))^2 + w^2)', 'independent', 'x', 'dependent', 'y', ...
    'coefficients', {'w', 'c'});
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.StartPoint = [w_guess, c_guess]; % put in guesses

%% ============================================================================================== %%
%% FIT AND EXTRACT
[fit_result, gof] = fit(x_data, y_data, ft, opts); % fit model to data
% extract the fit parameters:
w = fit_result.w;
c = fit_result.c;

%% ============================================================================================== %%
%% FITTED CURVE
y_fit = w^2 ./ ((2 * (x_data-c)).^2 + w^2);

%% ============================================================================================== %%
%% CHI-SQUARED PER DEGREE OF FREEDOM
resids = y_data - y_fit; 

end

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################