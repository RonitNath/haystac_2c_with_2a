function gain_out = FUNC_cavcorrect_smart_gain_v1(f_GHz, pow, auto_scale, cav_props) 

%% DESCRIPTION
% This function calls the FUNC_cavcorrect_smart_gain_v1 *after* performing a correction for a cavity
% seen in reflection. The cavity center freuqency, kappa_measure, and kappa_loss (FWHM's) must all
% be supplied in the cav_props object
%
%% HISTORY
% - v1 created by Dan Palken on 10/24/19
%
%% REFERENCES
% [1]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Gain Lorentzian Fit > Removing
% Cavity Dependance from Gain Measurements"
%

%% ============================================================================================== %%
%% INITIALIZE
Settings_and_Defaults_v4;
fig_folder = cav_props.fig_folder; 

%% ============================================================================================== %%
%% GLOBAL
global chi_mm2_fcn

%% ============================================================================================== %%
%% USER INPUT
% booleans:
plot_chi_mm2 = 0;
plot_cav_correction = 0;
save_figs = 01;

%% ============================================================================================== %%
%% FIX DIMENSIONS
if size(f_GHz) ~= size(pow)
    f_GHz = f_GHz.';
end

%% ============================================================================================== %%
%% UN-CAVITIZE
% see ref. [1] for the simply math on how to do this
k_d_GHz = (cav_props.kappa_meas_kHz - cav_props.kappa_loss_kHz)*1E-6; % difference
k_s_GHz = (cav_props.kappa_meas_kHz + cav_props.kappa_loss_kHz)*1E-6; % sum
det_GHz = f_GHz - cav_props.f_cav_GHz;

chi_mm2_pow = chi_mm2_fcn(k_d_GHz, k_s_GHz, det_GHz);
adj_pow = pow ./ chi_mm2_pow; 

%% ============================================================================================== %%
%% PLOT CAVITY REFLECTION
if plot_chi_mm2
    x_min = min(f_GHz); 
    x_max = max(f_GHz); 
    y_min = 0; 
    y_max = 1.01; 

    fig = figure; 
    hold on;
    plot(f_GHz, chi_mm2_pow, 'Color', SD.myblue);
    xlabel('$f_\mathrm{ENA}$ [GHz]');
    ylabel('$|\chi_{mm}|^2$');
    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    title([cav_props.label, ': Cavity Reflection']);
    hold off;
    
    if save_figs
        saveas(fig, [fig_folder, f, cav_props.label, '-Measurement Cavity Reflection.jpg'], 'jpeg');
        disp('cavity reflection figure saved successfully');
    end
end

%% ============================================================================================== %%
%% PLOT CAVITY CORRECTION
if plot_cav_correction    
    x_min = min(f_GHz); 
    x_max = max(f_GHz); 
    y_min = min([pow; adj_pow], [], 'all');
    y_max_bare = max([pow; adj_pow], [], 'all');
    y_max = 1.05 * y_max_bare;
    
    chi_mm2_scaled_pow = chi_mm2_pow * y_max_bare;
    
    fig = figure; 
    hold on;
    p3 = plot(f_GHz, chi_mm2_scaled_pow, 'Color', SD.myblue);
    p1 = plot(f_GHz, pow, 'Color', SD.neoncarrot, 'LineWidth', 4);
    p2 = plot(f_GHz, adj_pow, 'Color', SD.black, 'LineWidth', 1.5);
    xlabel('$f_\mathrm{ENA}$ [GHz]');
    ylabel('power');
    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    title([cav_props.label, ': Cavity Power-Gain Adjustmnet']);
    mylegend([p1, p2, p3], {'pow', 'pow$/|\chi_{mm}|^2$', ...
        [num2str(y_max_bare, 2), '$\times|\chi_{mm}|^2$']}, 'Location', 'Best');
    
    hold off;
    if save_figs
        saveas(fig, [fig_folder, f, cav_props.label, '-Measurement Cavity-Corrected Gain.jpg'], ...
            'jpeg');
        disp('cavity-corrected gain figure saved successfully');
    end
end

%% ============================================================================================== %%
%% CALL SMART GAIN FUCNTION
% final step is to pass along the relevant variables to the smart gain function
gain_out = FUNC_smart_gain_v1(f_GHz, adj_pow, auto_scale); % use the now-adjusted power

end

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################