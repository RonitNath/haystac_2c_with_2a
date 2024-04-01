function gain_out = FUNC_smart_gain_v1(f_GHz, pow, auto_scale) 

%% DESCRIPTION
% This function assumes a lot and gives back a little. Everything is returned in one output obejct. 
% Here are the assumptions: 
% - The data is a power gain Lorentzian likely (but not necessarily) from an ENA. 
% - No asusmptiuon of units on the power. 
% - Guesses for initial values are estiamted from the data, so that the user does not have to spend
% time estimating them out side the call to this function. 
%
%% HISTORY
% - v1 created by Dan Palken on 10/7/19
%
%% REFERENCES
% [1]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Gain Lorentzian Fit"
%

%% ============================================================================================== %%
%% USER INPUT
plot_gain_fit = 0;

%% ============================================================================================== %%
%% ROTATE
% make the inputs row vectors if they are column vectors
if size(f_GHz, 2) == 1
    f_GHz = f_GHz.';
end
if size(pow, 2) == 1
    pow = pow.';
end

% %% ============================================================================================== %%
% %% SCALE
% scale = mean(pow); 
% pow = pow / scale; 

%% ============================================================================================== %%
%% DERIVED QUANTITIES
n_fs = length(f_GHz);

%% ============================================================================================== %%
%% CHECKS
if n_fs ~= length(pow)
    error('Gain fitting length mismatch');
end

%% ============================================================================================== %%
%% GUESSES

guess.f_c_GHz = mean(f_GHz);
guess.G_max = max(pow) / min(pow);

% guessing the FWHM BW is a bit trickier
cent_ind = round(n_fs/2);
pow_l = pow(1:cent_ind);
pow_r = pow(cent_ind+1:end);

half_max = max(pow) / 2;
[~, ind_HM_l] = min(abs(pow_l - half_max));
[~, ind_HM_r] = min(abs(pow_r - half_max));
ind_HM_r = ind_HM_r + cent_ind; % map back on to the original vectors
f_HM_l_GHz = f_GHz(ind_HM_l);
f_HM_r_GHz = f_GHz(ind_HM_r);
guess.w_GHz = f_HM_r_GHz - f_HM_l_GHz;

guess.s_guess = (pow(1) + pow(end))/2; % not a great guess for the scale factor, but a simple one. 


[fit_result, ~, ~, ~] = FUNC_gain_fit_v8(f_GHz, pow, plot_gain_fit, guess.G_max, guess.w_GHz, ...
    guess.f_c_GHz, guess.s_guess, auto_scale);


%% ============================================================================================== %%
%% OUTPUT
gain_out.G_max = fit_result.G_max;
gain_out.w_GHz = fit_result.w;
if auto_scale
    gain_out.s = fit_result.s; 
else
    gain_out.s = 1; 
end
gain_out.f_c_GHz = fit_result.cent;

% we want to output a function that just takes an input frequency and gives the fitted value there.
% plot with s = 1 so we're actually plotting the gain:
gain_out.output_lor_fct = @(freq_GHz) 1 + (gain_out.w_GHz^2 * (gain_out.G_max-1)) ...
    ./ (4*(freq_GHz - gain_out.f_c_GHz).^2 + gain_out.w_GHz^2); 

% we also want to output a function that takes in a single power gain value and a detuning from the
% center, and using just the FWHM bandwidth w, gives back the value of the Lorentzian at any other
% point. See ref. [1], 10/9/19. This is useful when we are using the same bandwidth and one measured
% gain to get a single-quadrature gain Lorentzian
% G is power gain at detuning D. Dq is the detuning to query:
gain_out.query_fct = @(D_MHz, G, Dq_MHz) ...
    1 + (4*D_MHz^2 + (gain_out.w_GHz*1E3)^2)./(4*Dq_MHz.^2 + (gain_out.w_GHz*1E3)^2) * (G-1);

end

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################