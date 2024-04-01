function [] = FUNC_disp_abg_v5(UNIV, OP_data, abg, OP_ind, det_gain_MHz, save_figs, ...
    fig_folder)
%% DESCRIPTION
% Solves for some paramters of interest that the abg measurement gives us access to and makes
% relevant pplots
%
%% HISTORY
% - v1 created by Dan Palken on 7 Oct 2019
% - v2 created by DP on 21 Oct 2019
% - v3 created by DP on 22 Oct 2019
% - v4 created by DP on 6 Nov 2019. Considers S_c as constant
% - v5 created by DP on 7 Nov 2019. Calculates G_s using primarily b and g measurements
%
%% REFERENCES
% [1]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Squeezing/Hot Rod Calibration >
% Reducing Frequency Dependence"

%% ============================================================================================== %%
%% INITIAL
Settings_and_Defaults_v4;

%% ============================================================================================== %%
%% CONSTANTS
n_measures = 3; % a, b, g
ind_a = 1; ind_b = 2; ind_g = 3;

%% ============================================================================================== %%
%% USER INPUT
% boolean(s):
auto_scale = 01; % turn on if the ENA traces do not have a normalization curve already divided out
main_plot = 01;

freqs_combine = 100; % combine groups of this many frequencies

% FLAG - make sure this is at the start of the analysis band
S_c_avg_start_MHz = 0.025; % start at or near the start of the analysis band. 
S_c_avg_end_kappas = 0.5; % put end in units of cavity kappas (FWHMs)


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

if length(abg.f_GHz_a) == length(abg.f_GHz_b) & length(abg.f_GHz_a) == length(abg.f_GHz_g)
    disp('RF frequency vector lengths for abg are identical, as expected');
else
    error('non-identical RF freuqency vector lengths for abg');
end

%% ============================================================================================== %%
%% DERIVED QUANTITIES
n_ENA_fs = length(abg.f_GHz_a);
f_cav_t_avg_GHz = (OP_data.f_cav_t1_GHz + OP_data.f_cav_t2_GHz) / 2;
f_cav_r_GHz = OP_data.f_cav_r_GHz;
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

% construct object needed for the cavity-adjusted smart gain function
cav_props.kappa_meas_kHz = (OP_data.cav_kappa_ext_t1_kHz + OP_data.cav_kappa_ext_t2_kHz) / 2;
cav_props.kappa_loss_kHz = (OP_data.cav_kappa_loss_t1_kHz + OP_data.cav_kappa_loss_t2_kHz) / 2;
% cavity frequency will be the same for the different a/b/g
cav_props.f_cav_GHz = f_cav_r_GHz; % use reflection mean
cav_props.fig_folder = fig_folder; % also store a folder to save figures to in this object

N_Hp_fDep_1q_qta = interp1(OP_data.f4NH_IF_Hz, OP_data.N_Hp_qta/2, f_IF_Hz, [], 'extrap');

S_c_avg_end_MHz = S_c_avg_end_kappas * cav_kappa_tot_avg_kHz*1E-3;
[~, ind_S_c_avg_start] = min(abs(f_IF_Hz*1E-6 - S_c_avg_start_MHz));
[~, ind_S_c_avg_end] = min(abs(f_IF_Hz*1E-6 - S_c_avg_end_MHz));
S_c_avg_range = ind_S_c_avg_start:ind_S_c_avg_end;


%% ============================================================================================== %%
%% ENA SPECTROSCOPY
% the ENA measurements show the JPA 2-quad gain, and the cavity may appear as a small feature off to
% one side
P_ENA_a = abg.I_a.^2 + abg.Q_a.^2;
P_ENA_b = abg.I_b.^2 + abg.Q_b.^2;
P_ENA_g = abg.I_g.^2 + abg.Q_g.^2;

% the ENA freuqency vectors will not be the same for a as for b and g, note order is a, b, g:
cav_props.label = 'a'; % for plot labeling
g_out(1) = FUNC_cavcorrect_smart_gain_v1(abg.f_GHz_a, P_ENA_a, auto_scale, cav_props);
cav_props.label = 'b';
g_out(2) = FUNC_cavcorrect_smart_gain_v1(abg.f_GHz_b, P_ENA_b, auto_scale, cav_props);
cav_props.label = 'g';
g_out(3) = FUNC_cavcorrect_smart_gain_v1(abg.f_GHz_g, P_ENA_g, auto_scale, cav_props);

% calculate ENA gains at the relevant frequencies from the inbuilt functionality
ENA_G_a = g_out(1).output_lor_fct(abg.f_GHz_a);
ENA_G_b = g_out(2).output_lor_fct(abg.f_GHz_b);
ENA_G_g = g_out(3).output_lor_fct(abg.f_GHz_g);

G_A_1q_pow = [abg.G_a, abg.G_b, abg.G_g];
G_A_fdep_1q_pow = NaN(n_measures, n_f_IF); % order is a, b, g
for i = 1:n_measures
    G_A_fdep_1q_pow(i, :) = g_out(i).query_fct(det_gain_MHz, G_A_1q_pow(i), f_IF_Hz*1E-6);
end

%% ============================================================================================== %%
%% CALCULATE PARAMTERS OF INTEREST
% cav measurement port power reflection at the window where we perform our hot
% rod/squeezing calibration measurements:
chi_mm2_Wdet_pow = ((cav_kappa_diff_avg_kHz*1E3)^2 + 4*f_IF_Hz.^2) ./ ...
    ((cav_kappa_tot_avg_kHz*1E3)^2 + 4*f_IF_Hz.^2);

S_c_1q_qta = ((((N_Hp_fDep_1q_qta + G_A_fdep_1q_pow(ind_a,:) .* S_f_1q_qta) .* P_rat_ba - ...
    N_Hp_fDep_1q_qta) ./ G_A_fdep_1q_pow(ind_b,:)  - (1 - UNIV.cav_AMP_eta) * S_f_1q_qta) / ...
    UNIV.cav_AMP_eta - chi_mm2_Wdet_pow .* S_f_1q_qta) ./ (1 - chi_mm2_Wdet_pow);
% it is reasonable on physical grounds to assume that S_c is spectrally flat. Moreover S_c is in
% practice very hard to measure past the cavity linewidth, where the "signal" from the hot rod can
% be easily swamped by any systematic changes not accounted for in our simple model. Thus we use the
% average value of S_c over a range where it is both reasonable and relevant
S_c_avg_1q_qta = mean(S_c_1q_qta(S_c_avg_range));

G_s_pow = (((((N_Hp_fDep_1q_qta + G_A_fdep_1q_pow(ind_a,:) .* S_f_1q_qta) .* P_rat_ga - ...
    N_Hp_fDep_1q_qta) ./ G_A_fdep_1q_pow(ind_g,:) - (1 - UNIV.cav_AMP_eta) * S_f_1q_qta) / ...
    UNIV.cav_AMP_eta - (1 - chi_mm2_Wdet_pow) * S_c_avg_1q_qta) ./ chi_mm2_Wdet_pow - ...
    (1 - UNIV.SQ_cav_eta) * S_f_1q_qta) ./ (UNIV.SQ_cav_eta * S_f_1q_qta);

% include the full frequency dependence for comparison - this is not what we will use
G_s_full_fdep_pow = (((((N_Hp_fDep_1q_qta + G_A_fdep_1q_pow(ind_a,:) .* S_f_1q_qta) .* P_rat_ga - ...
    N_Hp_fDep_1q_qta) ./ G_A_fdep_1q_pow(ind_g,:) - (1 - UNIV.cav_AMP_eta) * S_f_1q_qta) / ...
    UNIV.cav_AMP_eta - (1 - chi_mm2_Wdet_pow) .* S_c_1q_qta) ./ chi_mm2_Wdet_pow - ...
    (1 - UNIV.SQ_cav_eta) * S_f_1q_qta) ./ (UNIV.SQ_cav_eta * S_f_1q_qta);

% there is an alternate, perhaps better (because it makes less use of alpha, the off-resonance 
% measruement) way to calculate G_s. See ref. [1]
z = N_Hp_fDep_1q_qta + G_A_fdep_1q_pow(ind_b,:) .* ((1-UNIV.cav_AMP_eta) * S_f_1q_qta + ...
    UNIV.cav_AMP_eta .* ((1-chi_mm2_Wdet_pow)*S_c_avg_1q_qta + chi_mm2_Wdet_pow * S_f_1q_qta));
G_s_alt_pow = ((((z .* P_rat_gb - N_Hp_fDep_1q_qta)./G_A_fdep_1q_pow(ind_g,:) - ...
    (1-UNIV.cav_AMP_eta)*S_f_1q_qta) / UNIV.cav_AMP_eta - (1 - chi_mm2_Wdet_pow) * S_c_avg_1q_qta) ./ ...
    chi_mm2_Wdet_pow - (1 - UNIV.SQ_cav_eta)*S_f_1q_qta) / (UNIV.SQ_cav_eta * S_f_1q_qta);

z_full_fdep = N_Hp_fDep_1q_qta + G_A_fdep_1q_pow(ind_b,:) .* ((1-UNIV.cav_AMP_eta) * S_f_1q_qta + ...
    UNIV.cav_AMP_eta .* ((1-chi_mm2_Wdet_pow).*S_c_1q_qta + chi_mm2_Wdet_pow * S_f_1q_qta));
G_s_alt_full_fdep_pow = ((((z .* P_rat_gb - N_Hp_fDep_1q_qta)./G_A_fdep_1q_pow(ind_g,:) - ...
    (1-UNIV.cav_AMP_eta)*S_f_1q_qta) / UNIV.cav_AMP_eta - (1 - chi_mm2_Wdet_pow) .* S_c_1q_qta) ./ ...
    chi_mm2_Wdet_pow - (1 - UNIV.SQ_cav_eta)*S_f_1q_qta) / (UNIV.SQ_cav_eta * S_f_1q_qta);

rod_partic = (S_c_avg_1q_qta - S_f_1q_qta) / (S_r_1q_qta - S_f_1q_qta);

%% ============================================================================================== %%
%% MAIN PLOT
if main_plot
    fig = figure('Position', [SD.default_left, SD.default_bot, SD.default_width, ...
        SD.default_height]);
    pos_subp1 = [0.1430, 0.40, 0.75, 0.53];
    subp1 = subplot('Position', pos_subp1);
    x_min = 0; % MHz
    x_max = f_plot_max_MHz; % MHz
    [~, ind_max] = min(abs(f_IF_Hz - f_plot_max_MHz*1E6));
    range = 1:ind_max;
    y_min = 0;
    P_abg_norm = P_abg ./ G_A_fdep_1q_pow;
    y_max = 1.05 * max(P_abg_norm(:, range), [], 'all');
    hold on;
    p1 = plot(f_IF_Hz*1E-6, P_abg_norm(ind_a,:), 'Color', SD.myred);
    p2 = plot(f_IF_Hz*1E-6, P_abg_norm(ind_b,:), 'Color', SD.myblue);
    p3 = plot(f_IF_Hz*1E-6, P_abg_norm(ind_g,:), 'Color', SD.mygreen);
    p4 = plot(S_c_avg_start_MHz*[1,1], [y_min, y_max], '--', 'Color', SD.mygray, ...
        'LineWidth', 1);
    
    ylabel('norm. PSD [a.u.]');
    title(['$\alpha\beta\gamma$ calibration measurements, OP ', num2str(OP_ind)]);
    mylegend([p1, p2, p3], {'$\alpha$ (SQ 0, res 0)', '$\beta$ (SQ 0, res 1)', ...
        '$\gamma$ (SQ 1, res 1)'}, 'Location', 'Best');
    xticks([]); % no ticks on x-axis
    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    
    subp2 = subplot('Position', [pos_subp1(1),  pos_subp1(2) - 0.3, pos_subp1(3), 0.3]);
    hold on;
    y_min = -0.5;
    y_max = 1; 
    % y_max = 1.1 * max([S_c_1q_qta(range), G_s_pow(range), rod_partic, S_r_1q_qta]);
    
    p1 = plot([x_min, x_max], S_r_1q_qta*[1, 1], '--', 'Color', SD.myorange, 'LineWidth', 1);
    p2 = plot([x_min, x_max], S_f_1q_qta*[1, 1], '--', 'Color', SD.myorange, 'LineWidth', 1);
    p3 = plot(f_IF_Hz*1E-6, S_c_1q_qta, 'Color', (SD.myorange + SD.white)/2, 'LineWidth', 2);
    p7 = plot(f_IF_Hz(S_c_avg_range)*1E-6, S_c_1q_qta(S_c_avg_range), ...
        'Color', SD.myorange, 'LineWidth', 2);
    p8 = plot([x_min, x_max], S_c_avg_1q_qta*[1,1], '-', 'Color', SD.myorange, 'LineWidth', 0.5);
    p9 = plot(f_IF_Hz*1E-6, G_s_full_fdep_pow, 'Color', SD.mypurple, 'LineWidth', 0.5);
    p10 = plot(f_IF_Hz*1E-6, G_s_alt_pow, 'Color', (SD.mypurple + 3*SD.black)/4, 'LineWidth', 2);
    p11 = plot(f_IF_Hz*1E-6, G_s_alt_full_fdep_pow, 'Color', (SD.mypurple + 3*SD.black)/4, ...
        'LineWidth', 0.5);
    
    p4 = plot(f_IF_Hz*1E-6, G_s_pow, 'Color', SD.mypurple, 'LineWidth', 2);
    p5 = plot([x_min, x_max], rod_partic*[1,1], 'Color', SD.goldenrod, 'LineWidth', 1);
    p6 = plot(S_c_avg_start_MHz*[1,1], [y_min, y_max], '--', 'Color', SD.mygray, ...
        'LineWidth', 1);
    p12 = plot([x_min, x_max], [0,0], 'Color', SD.mygray, 'LineWidth', 1);

    
    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    xlabel('$f_\mathrm{IF}$ [MHz]');
    ylabel('parameter');
%     mylegend([p3, p7, p8, p4, p10, p9, p11, p5, p2], {'$S_c$ (all $f$)', '$S_c$ (avg)', '$S_c$ (used)', ...
%         '$G_s$ (used)', '$G_s$ (alt)', '$G_s$ (all $f$)', '$G_s$ (alt all $f$)', '$p_r$', '$S_f,S_r$'}, 'Location', 'Best');
    hold off;
    if save_figs
        saveas(fig, [fig_folder, f, 'abg Calibrations.jpg'], 'jpeg');
        disp('abg calibrations figure saved successfully');
    end
end

end

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################