function [fit_result, y_fit, gof, chi2_per_DOF] = ...
    FUNC_centered_Lorentzian_fit_v1(x, y, w_guess, y_max_guess, weights)
%% DESCRIPTION: 
% Fits a Lorentzian known to be centered at 0. Fits for an overall scale factor y_max and a FWHM w
%
%% HISTORY:
% - v1 created on 4/9/18 by Dan Palken. Based closely upon FUNC_gain_fit_v5, but pared down 
% 
%% VARIABLE INFO:
% INPUT
% - x: independant variable data. Uncertainties assumed appropriately low
% - y: dependant variable data
% - y_max_guess: initial guess for the max height of the Lorentzian,
% - w_guess: inintial guess for the fit parameter w, corresponding to the full width at half max 
% (FWHM) of the Lorentzian
% - weights (corresponding to inverse varaiances of the power gains) to be used in the fit
% OUTPUT
% - fit_result: a fit object representing the fit
% - y_fit: fitted values at all the input locations
% - gof: structure with goodness-of fit info
% - chi2_per_DOF: chi-squared per degree of freedom for the fit
%

%% ============================================================================================== %%
%% USER INPUT
num_fit_params = 2; % make sure to update as program is altered

%% ============================================================================================== %%
%% PREPARE FIT
% assumes that the center of the Lorentzian is at 0
[x_data, y_data] = prepareCurveData(x, y);
ft = fittype('(y_max * w^2) / ((2 * x)^2 + w^2)', 'independent', 'x', 'dependent', 'y', ...
    'coefficients', {'w', 'y_max'});
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
opts.StartPoint = [w_guess, y_max_guess]; % put in guesses
opts.Weights = weights; % include weights

%% ============================================================================================== %%
%% FIT AND EXTRACT
[fit_result, gof] = fit(x_data, y_data, ft, opts); % fit model to data
% extract the fit parameters:
y_max = fit_result.y_max;
w = fit_result.w;

%% ============================================================================================== %%
%% FITTED CURVE
y_fit = (y_max * w^2) ./ ((2 * x_data).^2 + w^2);

%% ============================================================================================== %%
%% CHI-SQUARED PER DEGREE OF FREEDOM
resids = y_data - y_fit; 
n_DOFs = length(x) - num_fit_params; 
chi2 = sum((resids .^ 2) .* weights); 
chi2_per_DOF = chi2 / n_DOFs;
end

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################