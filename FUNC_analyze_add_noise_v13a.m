function [f_IF_Hz, N_H_1q_qta, G_H_1q_pow] = FUNC_analyze_add_noise_v13a(OP_data, OP_ind, ...
    T_VTS_K, VF, S_c_known_1q_qta, det_gain_MHz, det_norm_tone_MHz, det_dc2cav_MHz, T_f_K, ...
    lambda, rho, alpha, f_IF_Hz, S_out_norm_1q, f_ENA_cold_GHz, I_cold, Q_cold, f_ENA_hot_GHz, ...
    I_hot, Q_hot, save_figs, fig_folder, AMP_gain_meas)

%% DESCRIPTION
% handles the data from an added noise measurement. VTS temperature and and single-quadrature gain
% must have the same dimensionality (they are both independent variables, see ref. [1]), while the
% fridge temerature is assumed constant. The VTS-AMP losses are also assumed to be constant, even
% though **techncially** we could get it from this model with our two independent varaiables. Given
% how we take the measurement, that would be a bad idea
%
%% HISTORY
% - v1 created by Dan Palken on 9/23/19
% - v2 created by Dan Palken on 10/1/19
% - v3 created by Dan Palken on 10/4/19
% - v4 created by Dan Palken on 10/18/19
% - v5 created by Dan Palken on 10/21/19
% - v6 created by Dan Palken on 10/22/19
% - v7 created by Dan Palken on 11/11/19. Optionally removes gain-normalization tone. Emperically
% the several nearest neighbor points have lower power. Replace them all with the average of nearby
% points
% - v8 created by Dan Palken on 11/22/19. Allows for a custom S_c to be input at a known cavity
% detuning
% - v9 created by Dan Palken on 12/2/19. Uses smart cavity correction code to get gains
% - v10 created by Dan Palken on 12/2/19. Begins work on the value-function approach to optimizing
% calibration parameters
% - v11 created by DP on 12/4/19. Switches value function to taking in over gain ratio
% - v12 created by DP on 12/6/19
% - v13 created by DP on 12/31/19. Pairs with wkspc_construction_v14
% - v13a created by DP on 1/30/20. Pairs with wkspc_construction_v15a. Special changes flagged with 
% "CHG15a"
%
%% REFERENCES
% [1]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Y-Factor Measurement"
% [2]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Finitely Detuned Cavity &
% Calibration Measurements"
%

%% ============================================================================================== %%
%% INITIAL
Settings_and_Defaults_v4;
f = filesep;

%% ============================================================================================== %%
%% GLOBAL
global chi_mm2_fcn

%% ============================================================================================== %%
%% USER INPUT
% booleans:
plots_always_on = 0; % overrides functionality that turns plots off if no hot rod
plot_results = 0;
plot_ENA = 0;
plot_spectra = 0;
auto_scale = 01; % set true if 1 on the thing being fitted may not mean unit gain
remove_norm_tone = 01; % replace normalization tone with average of nearest neighbor points

% for averaging over the probe tone
freqs_dist = 3;
freqs_group = 5;

freqs_combine = 300; % combine groups of this many frequencies

%% ============================================================================================== %%
%% CHECKS
if lambda < 0 || lambda > 1
    error(['lambda = ', num2str(eta), ' is an unphysical value']);
end
if rho < 0 || rho > 1
    error(['rho = ', num2str(eta), ' is an unphysical value']);
end
if alpha < 0 || alpha > 1
    error(['alpha = ', num2str(eta), ' is an unphysical value']);
end
if T_f_K < 0
    error(['T_f_mK = ', num2str(T_f_K), ' is an unphysical value']);
end
if min(T_VTS_K) < 0
    error(['min(T_VTS_mK) = ', num2str(min(T_VTS_K)), ' is an unphysical value']);
end
if length(T_VTS_K) ~= size(VF.GA_addnoise_test_1q_pow, 2) % first non-1 dimensions should align 
    error('num-temps dimensions disagreement: T_VTS vs. G_A');
end
if f_ENA_cold_GHz ~= f_ENA_hot_GHz
    error('ENA frequency disagreement');
end
f_ENA_GHz = f_ENA_cold_GHz;
clear(varname(f_ENA_cold_GHz), varname(f_ENA_hot_GHz));

%% ============================================================================================== %%
%% DERIVED QUANTITIES
lambda_c = 1-lambda; % complement
rho_c = 1-rho; % complement
alpha_c = 1-alpha; % complement
f_meanRF_GHz = mean(f_ENA_GHz);
n_indep_OPs = size(S_out_norm_1q, 2); % number of points for the independent variables

% convert temperatures to single-quadrature spectral densities
S_VTS_1q_qta = FUNC_singlequad_S_of_T_f_v1(T_VTS_K, f_meanRF_GHz);
S_f_1q_qta = FUNC_singlequad_S_of_T_f_v1(T_f_K, f_meanRF_GHz);

G_A_test_1q_pow = VF.GA_addnoise_test_1q_pow;

n_GAs_test = size(G_A_test_1q_pow, 1); % how many gains we will try our value-function approach on

%% ============================================================================================== %%
%% SPECIAL FUNCTIONALITY
% the way this function will be called is such that on the first call S_c_known_1q_qta will be 0.
% When this is the case, it is because we are just in the context of looking for the first cavity
% HWHM of added noises. For that, we can reasonably safely assume we are far enough from the cavity
% that we can just say the fridge spectral density is basically what is coming out of the cavity
% (ref. [2]). The rest of the math will work itself out Include a message either way
if max(S_c_known_1q_qta,[],'all') == 0 && plots_always_on == 0
    disp('special value of S_c = 0 passed for added noise far from cavity: setting S_c = S_f');
    S_c_known_1q_qta = S_f_1q_qta;
    disp([tab, 'turning plots off']);
    plot_results = 0;
    plot_ENA = 0;
    plot_spectra = 0;
else
    disp('S_c not = 0 or plot overide on: N_H will be valid over full freqeuncy-dependence');
end

%% ============================================================================================== %%
%% REMOVE TONE
% remove the normalziation tone
[~, ind_norm] = min(abs(f_IF_Hz - det_norm_tone_MHz*1E6));
range_norm = (ind_norm-freqs_dist):(ind_norm+freqs_dist);
range_l = (ind_norm-freqs_dist-freqs_group):(ind_norm-freqs_dist-1);
range_r = (ind_norm+freqs_dist+1):(ind_norm+freqs_dist+freqs_group);
S_out_norm_1q(range_norm,:,:) = (ones(2*freqs_dist+1, 2, VF.n_hot_mults) .* ...
    mean([S_out_norm_1q(range_l,:,:); S_out_norm_1q(range_r,:,:)], 1));

%% ============================================================================================== %%
%% COARSEN
% the spectral resolution may be quite fine initially. We can combine the data in horizontal blocks
% to have fewer, more statistically reliable outputs
% first, kill the 1st element of the frequency/output vector, as it contibutes nothing since we
% subtract the mean is getting the PSD
f_IF_Hz = f_IF_Hz(2:end);
S_out_norm_1q = S_out_norm_1q(2:end,:,:);

n_f_IF = length(f_IF_Hz);
freqs_remain =  floor(n_f_IF / freqs_combine);
n_f_IF = freqs_remain * freqs_combine; % update number of IF frequencies

% cut end off to fascilitate reshaping:
f_IF_Hz = f_IF_Hz(1:n_f_IF);
S_out_norm_1q = S_out_norm_1q(1:n_f_IF,:,:);

% store a version before reshaping:
f_IF_preresahpe_Hz = f_IF_Hz;
S_out_norm_prereshape_1q = S_out_norm_1q;

shape = [freqs_combine, freqs_remain];
shape_out = [shape, n_indep_OPs, VF.n_hot_mults];
f_IF_HZ_reshp = reshape(f_IF_Hz, shape);
S_out_Gnorm_1q_reshp = NaN(shape_out); % init NaN
for i = 1:n_indep_OPs
    for j = 1:VF.n_hot_mults
        S_out_Gnorm_1q_reshp(:, :, i, j) = reshape(S_out_norm_1q(:, i, j), shape);
    end
end
f_IF_Hz = mean(f_IF_HZ_reshp);
S_out_norm_1q = permute(squeeze(mean(S_out_Gnorm_1q_reshp)),[2,1,3]);
n_f_IF = length(f_IF_Hz); % update number of IF frequencies once more

%% ============================================================================================== %%
%% ENA SPECTROSCOPY
% the ENA measurements show the JPA 2-quad gain, and the cavity may appear as a small feature off to
% one side
% CHG15a - lines commented out:
% P_tot_cold = I_cold.^2 + Q_cold.^2;
% P_tot_hot = I_hot.^2 + Q_hot.^2;

% construct object needed for the cavity-adjusted smart gain function:
cav_props.kappa_meas_kHz = (OP_data.cav_kappa_ext_t1_kHz + OP_data.cav_kappa_ext_t2_kHz) / 2;
cav_props.kappa_loss_kHz = (OP_data.cav_kappa_loss_t1_kHz + OP_data.cav_kappa_loss_t2_kHz) / 2;
% cavity frequency will be the same for the different a/b/g
cav_props.f_cav_GHz = (OP_data.f_cav_r_GHz + OP_data.f_cav_r2_GHz)/2; % use reflection mean
cav_props.fig_folder = fig_folder; % also store a folder to save figures to in this object

% CHG15a - we do not have ENA cals here, so use the one for the run (as with abgs): 
cav_props.label = 'run'; 
g_out_run = FUNC_cavcorrect_smart_gain_v1(AMP_gain_meas.AMPon_f_GHz, ...
    AMP_gain_meas.AMPon_meas_pow_big, auto_scale, cav_props);
g_out(1) = g_out_run; % set cold equal to run 
g_out(1).f_c_GHz = f_meanRF_GHz;
g_out(1).output_lor_fct = @(freq_GHz) 1 + (g_out(1).w_GHz^2 * (g_out(1).G_max-1)) ...
    ./ (4*(freq_GHz - g_out(1).f_c_GHz).^2 + g_out(1).w_GHz^2); 
g_out(1).query_fct = @(D_MHz, G, Dq_MHz) ...
    1 + (4*D_MHz^2 + (g_out(1).w_GHz*1E3)^2)./(4*Dq_MHz.^2 + (g_out(1).w_GHz*1E3)^2) * (G-1);
g_out(2) = g_out(1); % set hot equal to cold now

% calculate gains at the relevant frequencies from the inbuilt functionality:
G_cold = g_out(1).output_lor_fct(f_ENA_GHz);
G_hot = g_out(2).output_lor_fct(f_ENA_GHz);

%% ============================================================================================== %%
%% FREQUENCY-DEPENDENT GAIN
% next is to make the G_As_test_1q_pow object 4D. We want frequency-dependent values from gain_fct
% in the gain output object that will come from a gain fit of the ENA stuff. These will then apply
% at each frequency

G_A_test_fdep_1q_pow = NaN(n_GAs_test, n_indep_OPs, n_f_IF, VF.n_hot_mults);
G_A_test_normtonedet_1q_pow = NaN(n_GAs_test, n_indep_OPs, VF.n_hot_mults);
% calculate frequency dependent gains assuming cavity is well centered on JPA gain peak
for i = 1:n_indep_OPs
    for j = 1:n_GAs_test
        for k = 1:VF.n_hot_mults
            G_A_test_fdep_1q_pow(j, i, :, k) = g_out(i).query_fct(det_gain_MHz, ...
                G_A_test_1q_pow(j, i, k), f_IF_Hz*1E-6);
            
            % find gain at the detuning of the tone used to normalize - see below
            G_A_test_normtonedet_1q_pow(j, i, k) = g_out(i).query_fct(det_gain_MHz, ...
                G_A_test_1q_pow(j, i, k), det_norm_tone_MHz);
        end
    end
end

% Note that the output spectral densities have already been normalized by the gain at some given
% detuning. This is technically not the same as normalzing by the frequency dependent gains. We
% want to find the vlaues to correct by. The ones basically right at DC (i.e. closer than the
% detuning of the tone used to normalize them) will be > 1. Note that the detuning of the tone used
% for the single-quadrature gain measurement is *not* necessarily the same as that used to normalize
% all the spectra. For example, we may be talking about 400 kHz vs. 10 kHz detuning
rel_G_A_test = NaN(size(G_A_test_fdep_1q_pow));
for k = 1:VF.n_hot_mults
    rel_G_A_test(:,:,:,k) = G_A_test_fdep_1q_pow(:,:,:,k) ./ G_A_test_normtonedet_1q_pow(:,:,k);
end
S_out_test_norm_1q = NaN(size(rel_G_A_test));
for j = 1:n_GAs_test
    S_out_test_norm_1q(j,:,:,:) = S_out_norm_1q ./ squeeze(rel_G_A_test(j,:,:,:)); % renormalize
end

%% ============================================================================================== %%
%% GUESSES
guess_N_H_1q = 10; % 2X raw HEMT measurements from Yale in 2018

%% ============================================================================================== %%
%% NORMALIZE
S_out_test_norm_1q = S_out_test_norm_1q ./ mean(S_out_test_norm_1q, 2); % by individual means

%% ============================================================================================== %%
%% CALCULATE
% see ref. [2]:
% l ~ lambda, r ~ rho, a ~ alpha, x ~ (|chi_mm|^2+1)/2, Sh ~ VTS hot spec den, Sf ~ fridge spec den,
% Sc ~ cavity spec den, GAh ~ AMP again hot, GAc ~ AMP gain cold,
% rS = output gain-normalized spec den ratio hot:cold
% note we are explicitely assuming that Sf is also the cold load temperature
N_H_fct = @(l, r, a, x, Sh, Sf, Sc, GAh, GAc, rS) ...
    (((((l*Sh + (1-l)*Sf)*r + (1-r)*Sf)*x + Sc.*(1-x))*a + (1-a)*Sf) - ...
    rS.*((Sf*x + Sc.*(1-x))*a + (1-a)*Sf)) ./ (rS./GAc - 1./GAh);
G_H_fct = @(a, Soutnormc, Sf, Sc, x, NH, GAc) ...
    Soutnormc ./ (((Sf*x + Sc.*(1-x))*a + (1-a)*Sf) + NH./GAc);

S_out_norm_c = squeeze(S_out_test_norm_1q(:,1,:,:));
S_out_norm_h = squeeze(S_out_test_norm_1q(:,2,:,:));
S_rat = S_out_norm_h ./ S_out_norm_c; % defined as hot:cold
GA_hot = squeeze(G_A_test_fdep_1q_pow(:,2,:,:));
GA_cold = squeeze(G_A_test_fdep_1q_pow(:,1,:,:));
det_cav_close_MHz = det_dc2cav_MHz - f_IF_Hz*1E-6;
det_cav_far_MHz = det_dc2cav_MHz + f_IF_Hz*1E-6;

k_m_MHz = 1E-3*(OP_data.cav_kappa_ext_t1_kHz + OP_data.cav_kappa_ext_t2_kHz) / 2; % measruement port
k_l_MHz = 1E-3*(OP_data.cav_kappa_loss_t1_kHz + OP_data.cav_kappa_loss_t2_kHz) / 2; % loss port
k_s_MHz = k_m_MHz + k_l_MHz; % sum
k_d_MHz = k_m_MHz - k_l_MHz; % difference

chi_mm2_close_pow = chi_mm2_fcn(k_d_MHz, k_s_MHz, det_cav_close_MHz);
chi_mm2_far_pow = chi_mm2_fcn(k_d_MHz, k_s_MHz, det_cav_far_MHz);
x = (chi_mm2_close_pow + chi_mm2_far_pow)/2;

% directly calculated values: assumes cold index is first, hot second:
S_c_known_1q_qta = permute(S_c_known_1q_qta, [1,3,2]); % adds singleton dimension if 2D
N_H_1q_qta = N_H_fct(lambda, rho, alpha, x, S_VTS_1q_qta(2), S_f_1q_qta, S_c_known_1q_qta, ...
    GA_hot, GA_cold, S_rat);
G_H_1q_pow = G_H_fct(alpha, S_out_norm_c, S_f_1q_qta, S_c_known_1q_qta, x, N_H_1q_qta, GA_cold);

%% ============================================================================================== %%
%% PLOT RESULTS
if plot_results
    % for plotting, need to select some indicies (or make potentially many plots). A helpful thing
    % to do will be to choose the first, middle, and last indicies for the run/cal gain ratios 
    % tested. Similarly, for simplicity we will just plot the middle index for the h/c gain ratios
    if n_GAs_test <= 3
        plot_inds = 1:n_GAs_test;
    else
        mid_ind = floor((n_GAs_test + 1)/2); % floor only applies if even
        plot_inds = [1, mid_ind, n_GAs_test];
    end
    n_plot_inds = length(plot_inds);
    hc_plot_ind = floor((VF.n_hot_mults + 1)/2);
    
    S_out_fct = @(l, r, a, SV, Sf, Sc, NH, GA, GH, x) ...
        (((((l*SV + (1-l)*Sf)*r + (1-r)*Sf)*x + Sc.*(1-x))*a + (1-a)*Sf) + NH./GA).*GH;
    n_T_plt = 251;
    
    T_VTS_K_plt = linspace(0*min(T_VTS_K), 1.1*max(T_VTS_K), n_T_plt).';
    S_VTS_1q_plt = FUNC_singlequad_S_of_T_f_v1(T_VTS_K_plt, f_meanRF_GHz);
    G_A_test_1q_plt = NaN(n_GAs_test, n_T_plt, n_f_IF, VF.n_hot_mults);
    S_out_test_norm_fit = NaN(size(G_A_test_1q_plt));
    for j = 1:n_GAs_test
        G_A_test_1q_plt(j,:,:,:) = interp1(S_VTS_1q_qta, squeeze(G_A_test_fdep_1q_pow(j,:,:,:)), ...
            S_VTS_1q_plt, [], 'extrap');
        
        if length(S_c_known_1q_qta) == 1
            S_c_plt = S_c_known_1q_qta;
        else
            S_c_plt = S_c_known_1q_qta(j,:,:); 
        end
        S_out_test_norm_fit(j,:,:,:) = S_out_fct(lambda, rho, alpha, S_VTS_1q_plt, S_f_1q_qta, ...
            S_c_plt, N_H_1q_qta(j,:,:), squeeze(G_A_test_1q_plt(j,:,:,:)), ...
            G_H_1q_pow(j,:,:), x);
    end
    
    % renormalize by fitted gains:
    G_H_4renorm_1q_pow = permute(G_H_1q_pow, [1,4,2,3]); % adds signleton dimension if 3D
    S_out_test_Gnorm_fit = S_out_test_norm_fit ./ G_H_4renorm_1q_pow;
    S_out_test_Gnorm_1q = S_out_test_norm_1q ./ G_H_4renorm_1q_pow;
    
    plot_max_f_MHz = 3;
    plot_max_f_Hz = plot_max_f_MHz*1E6;
    [~, ind_max_f] = min(abs(f_IF_Hz - plot_max_f_Hz));
    
    % individual VTS fits figure:
    plot_every = 2;
    f_inds = 1:plot_every:ind_max_f;
    col_lo_f = SD.myred;
    col_hi_f = SD.myblue;
    
    x_min = min(T_VTS_K_plt)*1E3;
    x_max = max(T_VTS_K_plt)*1E3;
    
    for j = 1:n_plot_inds
        ind2plot = plot_inds(j);
        gain2plot_1q_dB = VF.GA_addnoise_scale_test_1q_dB(ind2plot); 
        title_str = ['$G_{A,1q}^\mathrm{hc-scale} = ', num2str(gain2plot_1q_dB, 3),'$ dB'];
        
        fig = figure;
        hold on;
        lgd_ord = [];
        
        
        for i = f_inds
            frac_hi_f = (i-1)/(ind_max_f-1);
            frac_lo_f = 1-frac_hi_f;
            col = col_lo_f*frac_lo_f + col_hi_f*frac_hi_f;
            
            plot([x_min, x_max], 0.25*[1,1], '-.', 'LineWidth', 1, 'Color', SD.mygray);
            plot(T_VTS_K*1E3, S_out_test_Gnorm_1q(ind2plot,:,i,hc_plot_ind), '.', ...
                'MarkerSize', 32, 'Color', SD.black);
            px = plot(T_VTS_K*1E3, S_out_test_Gnorm_1q(ind2plot,:,i,hc_plot_ind), '.', ...
                'MarkerSize', 25, 'Color', col);
            
            lgd_ord = [lgd_ord, px]; %#ok<AGROW>
            plot(T_VTS_K_plt*1E3, S_out_test_Gnorm_fit(ind2plot,:,i, hc_plot_ind), '--', ...
                'LineWidth', 2, 'Color', col);
        end
        
        mylegend([lgd_ord(1), lgd_ord(end)], {'low $f_\mathrm{IF}$', 'high $f_\mathrm{IF}$'}, ...
            'Location', 'Best');
        xlabel('$T_\mathrm{VTS}$ [mK]');
        ylabel('$S_\mathrm{out}/(G_A(T_\mathrm{VTS}) \times G_{H^\prime})$ [1q-qta]');
        title(['OP ', num2str(OP_ind), ', ', title_str]);
        xlim([x_min, x_max]);
        hold off;
        if save_figs
            saveas(fig, [fig_folder, f, 'OP ', num2str(OP_ind), ' Individual VTS Fits, ', ...
                num2str(gain2plot_1q_dB, 3), ' dB Gain.jpg'], 'jpeg');
            disp('individual VTS fits figure saved successfully');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % added noise fit results figure:
        col = SD.myorange;
        fig = figure;
        hold on;
        plot(f_IF_Hz(1:ind_max_f)*1E-6, N_H_1q_qta(ind2plot, 1:ind_max_f, hc_plot_ind), '.', ...
            'MarkerSize', 20, 'Color', col);
        plot(f_IF_Hz(1:ind_max_f)*1E-6, N_H_1q_qta(ind2plot, 1:ind_max_f, hc_plot_ind), '-', ...
            'LineWidth', 1.25, 'Color', col);
        title(['OP ', num2str(OP_ind), ' HEMT+ Added Noise, ', title_str]);
        xlabel('$f_\mathrm{IF}$ [MHz]');
        ylabel('$N_{H^\prime,1q}$ [qta]');
        hold off;
        if save_figs
            saveas(fig, [fig_folder, f, 'OP ', num2str(OP_ind), ...
                ' HEMT+ Added Noise Fit Results, ', num2str(gain2plot_1q_dB, 3), ...
                ' dB Gain.jpg'], 'jpeg');
            disp('HEMT+ added noise fit figure saved successfully');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gain fit results figure:
        col = SD.mypurple;
        fig = figure;
        hold on;
        plot(f_IF_Hz(1:ind_max_f)*1E-6, G_H_1q_pow(ind2plot,1:ind_max_f, hc_plot_ind) / ...
            mean(G_H_1q_pow(ind2plot,1:ind_max_f, hc_plot_ind)), ...
            '.', 'MarkerSize', 20, 'Color', col);
        plot(f_IF_Hz(1:ind_max_f)*1E-6, G_H_1q_pow(ind2plot,1:ind_max_f, hc_plot_ind) / ...
            mean(G_H_1q_pow(ind2plot,1:ind_max_f, hc_plot_ind)), ...
            '-', 'LineWidth', 1.25, 'Color', col);
        title(['OP ', num2str(OP_ind), ' HEMT+ Gain, ', title_str]);
        xlabel('$f_\mathrm{IF}$ [MHz]');
        ylabel('$G_{H^\prime,1q}/\mu$ [power units]');
        hold off;
        if save_figs
            saveas(fig, [fig_folder, f, 'OP ', num2str(OP_ind), ' Gain Fit Results, ', ...
                num2str(gain2plot_1q_dB, 3), ...
                ' dB Gain.jpg'], 'jpeg');
            disp('gain fit figure saved successfully');
        end
    end
end

%% ============================================================================================== %%
%% PLOT ENA MEASUREMENTS
if plot_ENA
    fig = figure;
    hold on;
    plot(f_ENA_GHz, G_cold, 'Color', SD.myblue);
    plot(f_ENA_GHz, G_hot, 'Color', SD.myred);
    xlabel('$f_\mathrm{ENA}$ [GHz]');
    ylabel('Direct Gain');
    title(['OP ', num2str(OP_ind)]);
    mylegend('cold', 'hot');
    hold off;
    if save_figs
        saveas(fig, [fig_folder, f, 'OP ', num2str(OP_ind), ...
            ' Added Noise ENA Measurements.jpg'], 'jpeg');
        disp('added noise ENA measurements figure saved successfully');
    end
end

%% ============================================================================================== %%
%% PLOT SPECTRA
if plot_spectra
    f_IF_preresahpe_MHz = f_IF_preresahpe_Hz*1E-6;
    
    x_min = 0.015;
    x_max = plot_max_f_MHz;
    [~, f_ind_min] = min(abs(f_IF_preresahpe_MHz - x_min));
    [~, f_ind_max] = min(abs(f_IF_preresahpe_MHz - x_max));
    range = f_ind_min:f_ind_max;
    
    % assumes just 2 VTS temperatures - or rather plots just the first two
    fig = figure;
    hold on;
    plot(f_IF_preresahpe_MHz(range), S_out_norm_prereshape_1q(range, 1, hc_plot_ind), ...
        'Color', SD.blue);
    plot(f_IF_preresahpe_MHz(range), S_out_norm_prereshape_1q(range, 2, hc_plot_ind), ...
        'Color', SD.red);
    title('Y-Factor Spectra');
    xlabel('$f_\mathrm{IF}$ [MHz]');
    ylabel('$S_\mathrm{out}$ [a.u.]');
    title(['OP ', num2str(OP_ind)]);
    mylegend({['$T_\mathrm{VTS} = ', num2str(T_VTS_K(1)*1E3), '$ mK'], ['$T_\mathrm{VTS} = ', ...
        num2str(T_VTS_K(2)*1E3), '$ mK']});
    xlim([x_min, x_max]);
    hold off;
    if save_figs
        saveas(fig, [fig_folder, f, 'OP ', num2str(OP_ind), ' VTS Spectra.jpg'], 'jpeg');
        disp('VTS spectra figure saved successfully');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot a version where the two lines are scaled to have the same mean
    mean_1 = mean(S_out_norm_prereshape_1q(range, 1));
    mean_2 = mean(S_out_norm_prereshape_1q(range, 2));
    mean_rat = mean_2/mean_1;
    
    fig = figure;
    hold on;
    plot(f_IF_preresahpe_MHz(range), mean_rat * S_out_norm_prereshape_1q(range, 1, hc_plot_ind), ...
        'Color', SD.blue);
    plot(f_IF_preresahpe_MHz(range), S_out_norm_prereshape_1q(range, 2, hc_plot_ind), ...
        'Color', SD.red);
    title(['OP ', num2str(OP_ind), ' Y-Factor Spectra']);
    xlabel('$f_\mathrm{IF}$ [MHz]');
    ylabel('$S_\mathrm{out}$ (scaled) [a.u.]');
    mylegend({['$T_\mathrm{VTS} = ', num2str(T_VTS_K(1)*1E3), '$ mK'], ['$T_\mathrm{VTS} = ', ...
        num2str(T_VTS_K(2)*1E3), '$ mK']});
    xlim([x_min, x_max]);
    hold off;
    
    if save_figs
        saveas(fig, [fig_folder, f, 'OP ', num2str(OP_ind), ' Rescaled VTS Spectra.jpg'], 'jpeg');
        disp('Rescaled VTS spectra figure saved successfully');
    end
end

end

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################