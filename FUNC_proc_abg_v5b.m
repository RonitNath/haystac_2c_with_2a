function [abg_out] = FUNC_proc_abg_v5b(UNIV, OP_data, abg, OP_ind, det_gain_MHz, VF, save_figs, ...
    fig_folder, all_plots_off)
%% DESCRIPTION
% Solves for some paramters of interest that the abg measurement gives us access to and makes
% relevant plots
%
%% HISTORY
% - v1 created by Dan Palken on 11/10/19, based directly on FUNC_disp_abg_v5. Now includes an
% output object
% - v2 created by DP on 22 Nov 2019. Includes an option fo turning all plots off
% - v3 created by DP on 3 Dec 2019. Begins work on the value-function approach to optimizing
% calibration parameters
% - v4 created by DP on 4 Dec 2019. Switches value function to taking in over gain ratio
% - v5 created by DP on 31 Dec 2019. Pairs with wkspc_construction_v14
% - v5a created by DP on 30 Jan 2020. Pairs with wkspc_construction_v15a. Special changes flagged 
% with "CHG15a"
% - v5b created by DP on 16 Feb 2020. Pairs with wkspc_construction_v15b. Special changes flagged 
% with "CHG15b"
%
%% REFERENCES
% [1]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Squeezing/Hot Rod Calibration >
% Reducing Frequency Dependence"

%% ============================================================================================== %%
%% INITIAL
Settings_and_Defaults_v4;

%% ============================================================================================== %%
%% GLOBAL
global chi_mm2_fcn

%% ============================================================================================== %%
%% CONSTANTS
n_measures = 3; % a, b, g
ind_a = 1; ind_b = 2; ind_g = 3;

%% ============================================================================================== %%
%% USER INPUT
% booleans:
auto_scale = 01; % turn on if the ENA traces do not have a normalization curve already divided out
main_plot = 0;

freqs_combine = 300; % combine groups of this many frequencies

% FLAG - make sure this is at the start of the analysis band
S_c_avg_start_MHz = 0.045; % start at or near the start of the analysis band.
S_c_avg_end_MHz = 0.4;

f_plot_max_MHz = 2.05;

%% ============================================================================================== %%
%% CHECKS
% make sure that the 3 abg frequency vectros are the same:
if abg.f_Hz_a == abg.f_Hz_b & abg.f_Hz_a == abg.f_Hz_g  %#ok<*AND2>
    disp('IF frequency vectors for abg are identical, as expected');
else
    error('non-identical IF freuqency vectors for abg');
end
f_IF_Hz = abg.f_Hz_a; %only need one name at this point for these identical vectors

% CHG15b - lines commented out: 
% if length(abg.f_GHz_a) == length(abg.f_GHz_b) & length(abg.f_GHz_a) == length(abg.f_GHz_g)
%     disp('RF frequency vector lengths for abg are identical, as expected');
% else
%     error('non-identical RF freuqency vector lengths for abg');
% end

%% ============================================================================================== %%
%% DERIVED QUANTITIES
if all_plots_off
    main_plot = 0;
end

% n_ENA_fs = length(abg.f_GHz_a); % CHG15b - commented out
f_cav_t_avg_GHz = (OP_data.f_cav_t1_GHz + OP_data.f_cav_t2_GHz) / 2;
f_cav_r_GHz = (OP_data.f_cav_r_GHz + OP_data.f_cav_r2_GHz)/2; % use reflection mean
S_r_1q_qta = FUNC_singlequad_S_of_T_f_v1(UNIV.T_rod_K, f_cav_t_avg_GHz);

n_f_IF = length(f_IF_Hz);
freqs_remain =  floor(n_f_IF / freqs_combine);
n_f_IF = freqs_remain * freqs_combine; % update number of IF frequencies
% cut end off to fascilitate reshaping:
f_IF_Hz = f_IF_Hz(1:n_f_IF);
quad_amp = NaN(n_f_IF, n_measures);
quad_amp(:,ind_a) = abg.quad_amp_a(1:n_f_IF);
quad_amp(:,ind_b) = abg.quad_amp_b(1:n_f_IF);
quad_amp(:,ind_g) = abg.quad_amp_g(1:n_f_IF);

% store a version before reshaping:
f_IF_preresahpe_Hz = f_IF_Hz;
quad_amp_prereshape_1q = quad_amp;

shape = [freqs_combine, freqs_remain];
shape_out = [shape, n_measures];
f_IF_HZ_reshp = reshape(f_IF_Hz, shape);
quad_amp_reshp = NaN(shape_out); % init NaN
for i = 1:n_measures
    quad_amp_reshp(:, :, i) = reshape(quad_amp(:, i), shape);
end
f_IF_Hz = mean(f_IF_HZ_reshp);
P_abg = squeeze(mean(quad_amp_reshp)).';
clear(varname(quad_amp));
n_f_IF = length(f_IF_Hz); % update number of IF frequencies once more

S_f_1q_qta = FUNC_singlequad_S_of_T_f_v1(UNIV.T_f_K, f_cav_t_avg_GHz);
P_rat_ba = P_abg(ind_b, :) ./ P_abg(ind_a, :);
P_rat_ga = P_abg(ind_g, :) ./ P_abg(ind_a, :);
P_rat_gb = P_abg(ind_g, :) ./ P_abg(ind_b, :);

cav_kappa_diff_t1_KHz = OP_data.cav_kappa_ext_t1_kHz - OP_data.cav_kappa_loss_t1_kHz;
cav_kappa_diff_t2_KHz = OP_data.cav_kappa_ext_t2_kHz - OP_data.cav_kappa_loss_t2_kHz;
cav_kappa_diff_avg_kHz = (cav_kappa_diff_t1_KHz + cav_kappa_diff_t2_KHz) / 2;
cav_kappa_tot_t1_KHz = OP_data.cav_kappa_ext_t1_kHz + OP_data.cav_kappa_loss_t1_kHz;
cav_kappa_tot_t2_KHz = OP_data.cav_kappa_ext_t2_kHz + OP_data.cav_kappa_loss_t2_kHz;
cav_kappa_tot_avg_kHz = (cav_kappa_tot_t1_KHz + cav_kappa_tot_t2_KHz) / 2;

% construct object needed for the cavity-adjusted smart gain function:
cav_props.kappa_meas_kHz = (OP_data.cav_kappa_ext_t1_kHz + OP_data.cav_kappa_ext_t2_kHz) / 2;
cav_props.kappa_loss_kHz = (OP_data.cav_kappa_loss_t1_kHz + OP_data.cav_kappa_loss_t2_kHz) / 2;
% cavity frequency will be the same for the different a/b/g
cav_props.f_cav_GHz = f_cav_r_GHz; % use reflection mean
cav_props.fig_folder = fig_folder; % also store a folder to save figures to in this object

n_GAs_test = size(VF.G_abg_test_1q_pow, 1); % how many gains we will try our value-function approach on
N_Hp_fDep_1q_qta = NaN(n_GAs_test, n_f_IF, VF.n_hot_mults);
for j = 1:n_GAs_test
    for k = 1:VF.n_hot_mults
        N_Hp_fDep_1q_qta(j,:,k) = interp1(OP_data.f4NH_IF_Hz, VF.N_Hp_test_1q_qta(j,:,k), ...
            f_IF_Hz, [], 'extrap');
    end
end

[~, ind_S_c_avg_start] = min(abs(f_IF_Hz*1E-6 - S_c_avg_start_MHz));
[~, ind_S_c_avg_end] = min(abs(f_IF_Hz*1E-6 - S_c_avg_end_MHz));
S_c_avg_range = ind_S_c_avg_start:ind_S_c_avg_end;


%% ============================================================================================== %%
%% ENA SPECTROSCOPY
% the ENA measurements show the JPA 2-quad gain, and the cavity may appear as a small feature off to
% one side

% CHG15a - there are no gain traces for the abgs, so use the run gain for all three, but move the
% frequency for a
% order is abg:
cav_props.label = 'b';
g_out(2) = FUNC_cavcorrect_smart_gain_v1(abg.AMPon_f_GHz, abg.AMPon_meas_pow_big, auto_scale, cav_props);
g_out(3) = g_out(2);
g_out(1) = g_out(2);

% CHG15a - re-define the g_out(1) functions about the correct frequency 
% CHG15b - this line gives the proper center frequency for the 'a' calibration:
g_out(1).f_c_GHz = OP_data.f_cav_t2_GHz - UNIV.det_norm_tone_MHz*1E-3; 
g_out(1).output_lor_fct = @(freq_GHz) 1 + (g_out(1).w_GHz^2 * (g_out(1).G_max-1)) ...
    ./ (4*(freq_GHz - g_out(1).f_c_GHz).^2 + g_out(1).w_GHz^2); 
g_out(1).query_fct = @(D_MHz, G, Dq_MHz) ...
    1 + (4*D_MHz^2 + (g_out(1).w_GHz*1E3)^2)./(4*Dq_MHz.^2 + (g_out(1).w_GHz*1E3)^2) * (G-1);

% calculate ENA gains at the relevant frequencies from the inbuilt functionality
% CHG15b - lines commented out
% ENA_G_a = g_out(1).output_lor_fct(abg.f_GHz_a);
% ENA_G_b = g_out(2).output_lor_fct(abg.f_GHz_b);
% ENA_G_g = g_out(3).output_lor_fct(abg.f_GHz_g);

G_A_fdep_1q_pow = NaN(n_GAs_test, n_measures, n_f_IF); % order is a, b, g
for j = 1:n_GAs_test
    for i = 1:n_measures
        G_A_fdep_1q_pow(j,i,:) = g_out(i).query_fct(det_gain_MHz, VF.G_abg_test_1q_pow(j,i), ...
            f_IF_Hz*1E-6);
    end
end

%% ============================================================================================== %%
%% CALCULATE PARAMTERS OF INTEREST
% cav measurement port power reflection at the window where we perform our hot rod/squeezing
% calibration measurements:
chi_mm2_Wdet_pow = chi_mm2_fcn(cav_kappa_diff_avg_kHz, cav_kappa_tot_avg_kHz, f_IF_Hz*1E-3);

S_c_test_1q_qta = ((((N_Hp_fDep_1q_qta + squeeze(G_A_fdep_1q_pow(:,ind_a,:)) .* S_f_1q_qta) .* ...
    P_rat_ba - N_Hp_fDep_1q_qta) ./ squeeze(G_A_fdep_1q_pow(:,ind_b,:))  - ...
    (1 - UNIV.cav_AMP_eta) * S_f_1q_qta) / UNIV.cav_AMP_eta - chi_mm2_Wdet_pow .* S_f_1q_qta) ./ ...
    (1 - chi_mm2_Wdet_pow);
% it is reasonable on physical grounds to assume that S_c is spectrally flat. Moreover S_c is in
% practice very hard to measure past the cavity linewidth, where the "signal" from the hot rod can
% be easily swamped by any systematic changes not accounted for in our simple model. Thus we use the
% average value of S_c over a range where it is both reasonable and relevant
S_c_test_avg_1q_qta = mean(S_c_test_1q_qta(:,S_c_avg_range,:),2);

% Use the "alternate" (likely better) way to calculate G_s of ref. [1]
z_test = N_Hp_fDep_1q_qta + squeeze(G_A_fdep_1q_pow(:,ind_b,:)) .* ((1-UNIV.cav_AMP_eta) * ...
    S_f_1q_qta + UNIV.cav_AMP_eta .* (S_c_test_avg_1q_qta.*(1-chi_mm2_Wdet_pow) + ...
    chi_mm2_Wdet_pow * S_f_1q_qta));
G_s_test_pow = ((((z_test .* P_rat_gb - N_Hp_fDep_1q_qta)./squeeze(G_A_fdep_1q_pow(:,ind_g,:)) - ...
    (1-UNIV.cav_AMP_eta)*S_f_1q_qta) / UNIV.cav_AMP_eta - S_c_test_avg_1q_qta .* ...
    (1 - chi_mm2_Wdet_pow)) ./ chi_mm2_Wdet_pow - (1 - UNIV.SQ_cav_eta)*S_f_1q_qta) / ...
    (UNIV.SQ_cav_eta * S_f_1q_qta);

% this uses the fully frequency-dependent version of S_c:
z_test_full_fdep = N_Hp_fDep_1q_qta + squeeze(G_A_fdep_1q_pow(:,ind_b,:)) .* ...
    ((1-UNIV.cav_AMP_eta) * S_f_1q_qta + UNIV.cav_AMP_eta .* ...
    ((1-chi_mm2_Wdet_pow).* S_c_test_1q_qta + chi_mm2_Wdet_pow * S_f_1q_qta));
G_s_test_full_fdep_pow = ((((z_test .* P_rat_gb - N_Hp_fDep_1q_qta) ./ ...
    squeeze(G_A_fdep_1q_pow(:,ind_g,:)) - (1-UNIV.cav_AMP_eta)*S_f_1q_qta) / UNIV.cav_AMP_eta - ...
    (1 - chi_mm2_Wdet_pow) .* S_c_test_1q_qta) ./ chi_mm2_Wdet_pow - ...
    (1 - UNIV.SQ_cav_eta)*S_f_1q_qta) / (UNIV.SQ_cav_eta * S_f_1q_qta);

rod_partic_test = squeeze((S_c_test_avg_1q_qta - S_f_1q_qta) / (S_r_1q_qta - S_f_1q_qta));

%% ============================================================================================== %%
%% OUTPUT
S_c_test_avg_1q_qta = squeeze(S_c_test_avg_1q_qta);
% constrcut output object:
abg_out.S_c_test_1q_qta = S_c_test_avg_1q_qta;
abg_out.G_s_test_pow = G_s_test_pow;
abg_out.rod_partic_test = rod_partic_test;
abg_out.f_IF_Hz = f_IF_Hz;

%% ============================================================================================== %%
%% MAIN PLOT
if main_plot
    if n_GAs_test <= 3
        plot_inds = 1:n_GAs_test;
    else
        mid_ind = floor((n_GAs_test + 1)/2); % floor only applies if even
        plot_inds = [1, mid_ind, n_GAs_test];
    end
    n_plot_inds = length(plot_inds);
    hc_plot_ind = floor((VF.n_hot_mults + 1)/2);
    
    for j = 1:n_plot_inds
        ind2plot = plot_inds(j);
        gain2plot_1q_dB = VF.GA_addnoise_scale_test_1q_dB(ind2plot);
        title_str = ['$G_{A,1q}^\mathrm{hc-scale} = ', num2str(gain2plot_1q_dB, 3),'$ dB'];
        
        fig = figure('Position', [SD.default_left, SD.default_bot, SD.default_width, ...
            SD.default_height]);
        pos_subp1 = [0.1430, 0.40, 0.75, 0.53];
        subp1 = subplot('Position', pos_subp1);
        x_min = 0; % MHz
        x_max = f_plot_max_MHz; % MHz
        [~, ind_max] = min(abs(f_IF_Hz - f_plot_max_MHz*1E6));
        range = 1:ind_max;
        y_min = 0;
        P_abg_norm = P_abg ./ squeeze(G_A_fdep_1q_pow(ind2plot,:,:));
        y_max = 1.05 * max(P_abg_norm(:, range), [], 'all');
        hold on;
        p1 = plot(f_IF_Hz*1E-6, P_abg_norm(ind_a,:), 'Color', SD.myred);
        p2 = plot(f_IF_Hz*1E-6, P_abg_norm(ind_b,:), 'Color', SD.myblue);
        p3 = plot(f_IF_Hz*1E-6, P_abg_norm(ind_g,:), 'Color', SD.mygreen);
        p4 = plot(S_c_avg_start_MHz*[1,1], [y_min, y_max], '--', 'Color', SD.mygray, ...
            'LineWidth', 1);
        
        ylabel('norm. PSD [a.u.]');
        title(['$\alpha\beta\gamma$: OP ', num2str(OP_ind), ', ', title_str]);
        mylegend([p1, p2, p3], {'$\alpha$ (SQ 0, res 0)', '$\beta$ (SQ 0, res 1)', ...
            '$\gamma$ (SQ 1, res 1)'}, 'Location', 'Best');
        xticks([]); % no ticks on x-axis
        xlim([x_min, x_max]);
        ylim([y_min, y_max]);
        
        subp2 = subplot('Position', [pos_subp1(1),  pos_subp1(2) - 0.3, pos_subp1(3), 0.3]);
        hold on;
        y_min = -0.5;
        y_max = 1;
        
        p1 = plot([x_min, x_max], S_r_1q_qta*[1, 1], '--', 'Color', SD.myorange, 'LineWidth', 1);
        p2 = plot([x_min, x_max], S_f_1q_qta*[1, 1], '--', 'Color', SD.myorange, 'LineWidth', 1);
        p3 = plot(f_IF_Hz*1E-6, S_c_test_1q_qta(ind2plot,:,hc_plot_ind), ...
            'Color', (SD.myorange + SD.white)/2, 'LineWidth', 2);
        p7 = plot(f_IF_Hz(S_c_avg_range)*1E-6, ...
            S_c_test_1q_qta(ind2plot, S_c_avg_range, hc_plot_ind), 'Color', SD.myorange, ...
            'LineWidth', 2);
        p8 = plot([x_min, x_max], S_c_test_avg_1q_qta(ind2plot, hc_plot_ind)*[1,1], '-', ...
            'Color', SD.myorange, 'LineWidth', 0.5);
        p10 = plot(f_IF_Hz*1E-6, G_s_test_pow(ind2plot,:,hc_plot_ind), 'Color', SD.mypurple, ...
            'LineWidth', 2);
        p11 = plot(f_IF_Hz*1E-6, G_s_test_full_fdep_pow(ind2plot,:,hc_plot_ind), ...
            'Color', SD.mypurple, 'LineWidth', 0.5);
        
        p5 = plot([x_min, x_max], rod_partic_test(ind2plot)*[1,1], 'Color', SD.goldenrod, ...
            'LineWidth', 1);
        p6 = plot(S_c_avg_start_MHz*[1,1], [y_min, y_max], '--', 'Color', SD.mygray, ...
            'LineWidth', 1);
        p12 = plot([x_min, x_max], [0,0], 'Color', SD.mygray, 'LineWidth', 1);
        
        
        xlim([x_min, x_max]);
        ylim([y_min, y_max]);
        xlabel('$f_\mathrm{IF}$ [MHz]');
        ylabel('parameter');
        mylegend([p3, p7, p8, p10, p11, p5, p2], {'$S_c$ (all $f$)', '$S_c$ (avg)', ...
            '$S_c$ (used)', '$G_s$', '$G_s$ (full $f$-dep)', '$p_r$', '$S_f,S_r$'}, ...
            'Location', 'Best');
        hold off;
        if save_figs
            saveas(fig, [fig_folder, f, 'abg Calibrations, ', ...
                num2str(gain2plot_1q_dB, 3), ' dB Gain.jpg'], 'jpeg');
            disp('abg calibrations figure saved successfully');
        end
    end
end

end

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################