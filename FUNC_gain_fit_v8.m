function [fit_result, G_fit, gof, resids] = ...
    FUNC_gain_fit_v8(freq, G, plot_fit, G_max_guess, w_guess, cent_guess, s_guess, auto_scale)

% DESCRIPTION: 
% performs a Lorentzian fit of a frequency-dependent gain. Allows user to
% enter bounds and initial guesses for the fit parameters. 
%
%
% HISTORY:
% Program created as a modified version of Max's by Dan Palken on 1 Oct.
% 2017. First version is optimized for the case where the central frequency
% is known (it is the pump frequency). It is still a fit parameter, but the
% bounds are made narrow.
% v2 created 22 Jan 2017. Plotting boolean now an input option. 
% - y0 and xc fixed - i.e. no long fit parameters
% v3 created 24 Jan. 2017 by DP
% - xc no longer fixed - had forgotten about Kerr phenomena
% - likely v2 is the superior version to use, as I do not believe Kerr 
%   shifting applies to this physics. 
% v4 created 25 Jan. 2017 by DP
% - back in the style of v3 with fixed x_c
% - takes in guesses for a and w
% - if guesses are out of range, uses default guesses
% v5 created 6 Oct. 2017 by DP
% - recasts the fit in terms of more experimentally accessisble/easily
%   interpretable parameters.
%
%%
%%  Data for fit:
%%      - x-input: freq: frequency in GHz
%%      - y-input: G: gain in power units
%%
%%  Additonal input
%%      - plot_fit: boolean for plotting the data and fit 
%%      - G_max_guess: initial guess for the fit parameter G_max,
%%        corresponding to the maximum power gain: the ratio of the Lotentzian 
%%        at 0 vs infinite detuning. 
%%      - w_guess: inintial guess for the fit parameter w, corresponding to
%%        the full width at half max (FWHM) of the Lorentzian, in GHz
%%      - cent_guess: initial guess for the fit paramter cent, corresponding
%%        to the center frequency of the Lorentzian
%%      - s_guess: initial guess for the scale paramter s, corresponding to the 
%%        value of the Lorentzian at infinite detuning
%%      - auto_scale: boolean. If off, s = 1
%%
%%  Output:
%%      - fit_result: a fit object representing the fit.
%%      - G_fit: the fitted values for the gain at all the input 
%%        frequencies
%%      - gof: structure with goodness-of fit info.
%%      - resids: residuals from the fit at each input frequency, defined 
%%        as measured minus fitted gains
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE
Settings_and_Defaults_v4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DERIVED QUANTITIES
[xData, yData] = prepareCurveData(freq,G);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCALE
scale = mean(G); 
if auto_scale % only bother scaling down the numbers if we really dont know the scale
    G = G / scale; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP FIT
if auto_scale
    ft = fittype('s*(1 + (w^2 * (G_max-1)) / (4*(freq-cent)^2 + w^2))', 'independent', 'freq', ...
        'dependent', 'G');
else 
    ft = fittype('(1 + (w^2 * (G_max-1)) / (4*(freq-cent)^2 + w^2))', 'independent', 'freq', ...
        'dependent', 'G');
end

opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';

% bounds and guesses:
G_max_lower = 0.8; % power units
G_max_start = 1E2; % power units; default
G_max_upper = 1E6; % power units
w_lower = 0.0001; % GHz
w_start = 0.14; % GHz; default
w_upper = 1; % GHz
cent_lower = min(freq);
cent_start = median(freq);
cent_upper = max(freq);
s_lower = 0; 
s_start = 1;
s_upper = Inf;

% adjust default start values if input gueeses are in range:
if G_max_guess >= G_max_lower && G_max_guess <= G_max_upper
    G_max_start = G_max_guess; % power units
end
if cent_guess >= cent_lower && cent_guess <= cent_upper
    cent_start = cent_guess; % GHz
end
if w_guess >= w_lower && w_guess <= w_upper
    w_start = w_guess; % GHz
end
if s_guess >= s_lower && s_guess <= s_upper
    s_start = s_guess; 
end

if auto_scale
    opts.Lower = [G_max_lower, cent_lower, s_lower, w_lower];
    opts.StartPoint = [G_max_start, cent_start, s_start, w_start];
    opts.Upper = [G_max_upper, cent_upper, s_upper, w_upper];
else
    opts.Lower = [G_max_lower, cent_lower, w_lower];
    opts.StartPoint = [G_max_start, cent_start, w_start];
    opts.Upper = [G_max_upper, cent_upper, w_upper];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORM FIT
[fit_result, gof] = fit(xData, yData, ft, opts); % fit model to data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT
if plot_fit
    ws = warning('off', 'all');  

    % Plot fit with data.
    figure('Name', 'Lorentzian Gain Fit');
    plt = plot(fit_result, xData, yData);
    movegui('northwest');
    legend(plt, 'gain vs. freq', 'Lorentzian Gain Fit', 'Location', ...
        'NorthEast');
    xlabel("freq [GHz]");
    ylabel ("power gain");
    box on;
    ws = warning('on', 'all');  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FREE PARAMETERS 
G_max = fit_result.G_max;
w = fit_result.w;
cent = fit_result.cent;
if auto_scale
    s = fit_result.s * scale;
else
    s = 1; % no auto-scaling: same as s = 1
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FITTED CURVE
G_fit = s*(1 + (w^2 * (G_max-1)) ./ (4*(freq-cent).^2 + w^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESIDUALS 
G = G * scale; 
resids = G - G_fit;

end

% #########################################################################
% ############################ END OF FUNCTION ############################
% #########################################################################