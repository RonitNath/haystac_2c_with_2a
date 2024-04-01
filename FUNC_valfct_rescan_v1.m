function params_select = FUNC_valfct_rescan_v1(UNIV, OP_data, OP_ind, abg_test_out, VF, g_out, ...
    save_figs, fig_folder)
%% DESCRIPTION
% Takes in test number-dependent values for the HEMT added noise N_H, the cavity spectral density
% S_c, and the squeezing G_s, as well as the run gain G_A. Following the logic developed in refs.
% [1] and [2], evaluates the outputs against selected criteria and chooses that which minimizes the
% value function proportional to scan rate (to be conservative)
%
%% HISTORY
% - v1 created by Dan Palken on 24 Mar 2020. Based directly off of FUNC_valfct_v6
%
%% REFERENCES
% [1]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Relating the Gains"
% [2]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Relating the Gains > 
% Choosing the Optimal Hot/Cold Cal Gain Ratio"
%

%% ============================================================================================== %%
%% INITIAL
Settings_and_Defaults_v4;

%% ============================================================================================== %%
%% GLOBAL
global chi_mm2_fcn

%% ============================================================================================== %%
%% CONSTANTS
f = filesep;

%% ============================================================================================== %%
%% USER INPUT
% boolean(s):
plot_valfct = 01;
plot_params = 01;
plot_valid_param_space = 01; % plots a 2D map of the validity in parameter space

n_f_IFs = 250;
f_IF_min_int_kHz = 11;
f_IF_max_int_kHz = 1311;

%% ============================================================================================== %%
%% DERIVED QUANTITIES
fIF4NH_Hz = OP_data.f4NH_IF_Hz;
min_f4NH_Hz = min(fIF4NH_Hz);
max_f4NH_Hz = max(fIF4NH_Hz);
fIF4abg_Hz = abg_test_out.f_IF_Hz;
min_f4abg_Hz = min(fIF4abg_Hz);
max_f4abg_Hz = max(fIF4abg_Hz);
min_fIF_Hz = max(min_f4NH_Hz, min_f4abg_Hz); % min f to use is max of individual mins
max_fIF_Hz = min(max_f4NH_Hz, max_f4abg_Hz); % max f to use is min of individual maxs
f_IF_Hz = linspace(min_fIF_Hz, max_fIF_Hz, n_f_IFs); % f vector that we will use
[~, ind_f_min] = min(abs(f_IF_min_int_kHz*1E3 - f_IF_Hz));
[~, ind_f_max] = min(abs(f_IF_max_int_kHz*1E3 - f_IF_Hz));
range_ind_f = ind_f_min:ind_f_max;
f_IF_Hz = f_IF_Hz(range_ind_f); % narrow f down to this range
f_IF_min_int_kHz = min(f_IF_Hz)*1E-3; % re-adjust
f_IF_max_int_kHz = max(f_IF_Hz)*1E-3; % re-adjust
n_f_IFs = length(f_IF_Hz); % re-adjust 
n_GAs = VF.n_GAs;

f_cav_t_avg_GHz = (OP_data.f_cav_t1_GHz + OP_data.f_cav_t2_GHz) / 2;
S_f_1q_qta = FUNC_singlequad_S_of_T_f_v1(UNIV.T_f_K, f_cav_t_avg_GHz);

cav_kappa_diff_t1_KHz = OP_data.cav_kappa_ext_t1_kHz - OP_data.cav_kappa_loss_t1_kHz;
cav_kappa_diff_t2_KHz = OP_data.cav_kappa_ext_t2_kHz - OP_data.cav_kappa_loss_t2_kHz;
cav_kappa_diff_avg_kHz = (cav_kappa_diff_t1_KHz + cav_kappa_diff_t2_KHz) / 2;
cav_kappa_tot_t1_KHz = OP_data.cav_kappa_ext_t1_kHz + OP_data.cav_kappa_loss_t1_kHz;
cav_kappa_tot_t2_KHz = OP_data.cav_kappa_ext_t2_kHz + OP_data.cav_kappa_loss_t2_kHz;
cav_kappa_tot_avg_kHz = (cav_kappa_tot_t1_KHz + cav_kappa_tot_t2_KHz) / 2;

% assign concise variable names to the parameters which will show up in the nasty expression for
% the visibility. It is understood that spectral densities, added noises, and gains are single-quad
% and in linear-power units:
x = chi_mm2_fcn(cav_kappa_diff_avg_kHz, cav_kappa_tot_avg_kHz, f_IF_Hz*1E-3);
xc = 1-x;
y = 1 ./ ((cav_kappa_tot_avg_kHz^2)/4 + (f_IF_Hz*1E-3).^2); % |chi_ma|^2 with constants pulled out
Sf = FUNC_singlequad_S_of_T_f_v1(UNIV.T_f_K, f_cav_t_avg_GHz);
r = UNIV.SQ_cav_eta;
rc = 1-r;
a = UNIV.cav_AMP_eta;
ac = 1-a;
% GS will be all ones for the rescan version
GS = permute(interp1(fIF4abg_Hz, permute(abg_test_out.G_s_test_pow,[2,1,3]), f_IF_Hz),[2,1,3]);
NH = permute(interp1(fIF4NH_Hz, permute(VF.N_Hp_test_1q_qta,[2,1,3]), f_IF_Hz),[2,1,3]);
Sc = permute(abg_test_out.S_c_test_1q_qta,[1,3,2]); % adds singleton dimension if 1D
GA = g_out.query_fct(0, VF.GA_spec_fixed_1q_pow, f_IF_Hz*1E-6);

%% ============================================================================================== %%
%% EVALUATE
% calculate the visibility (or the thing proprtional to it) from ref. [1]
vis_num = y .* GA; % numerator
vis_den = (((GS*Sf*r + rc*Sf).*x + Sc.*xc)*a + ac*Sf).*GA + NH; % denominator
vis = vis_num ./ vis_den;
vis2 = vis.^2;
vis2int = vis2;

Sc = squeeze(Sc); % remove singleton dimension

R = squeeze(trapz(f_IF_Hz, vis2int, 2)); % value function goes as scan rate. Minimize ~ conservative

%% ============================================================================================== %%
%% VALIDITY
NHmin_crit = (NH >= VF.crit.NHmin_1q);
NHmin_crit_met = squeeze(sum(NHmin_crit, 2)/n_f_IFs >= VF.crit.fracmatch);
NH_crit_met = NHmin_crit_met; % no max criterion

GSmin_crit = (GS >= VF.crit.Gs_min_1q_pow); % this will be auto-satisfied for the rescan version
GSmin_crit_met = squeeze(sum(GSmin_crit, 2)/n_f_IFs >= VF.crit.fracmatch);
GSmax_crit = (GS <= VF.crit.Gs_max_1q_pow); % this will be auto-satisfied for the rescan version
GSmax_crit_met = squeeze(sum(GSmax_crit, 2)/n_f_IFs >= VF.crit.fracmatch);
GS_crit_met = (GSmin_crit_met & GSmax_crit_met);

Scmin_crit_met = (Sc >= VF.crit.Scmin_1q_qta);
Scmax_crit_met = (Sc <= VF.crit.Scmax_1q_qta);
Sc_crit_met = Scmin_crit_met & Scmax_crit_met;

all_crit_met = NH_crit_met & GS_crit_met & Sc_crit_met;

% now do a similar but separate check for if N_H, G_S, or S_c ever go zub-zero, as these are more
% physcially stringent conditions, and we observe some simularities e.g. for G_S <= 0
NH_sub0 = squeeze(logical(ceil(sum(NH <= 0, 2)/n_f_IFs)));
GS_sub0 = squeeze(logical(ceil(sum(GS <= 0, 2)/n_f_IFs)));
Sc_sub0 = (Sc <= 0);
any_sub0 = NH_sub0 | GS_sub0 | Sc_sub0;

% going below zero for our paramters is flatly dissalowed:
NH_crit_met = NH_crit_met & ~ NH_sub0;
GS_crit_met = GS_crit_met & ~ GS_sub0;
Sc_crit_met = Sc_crit_met & ~ Sc_sub0;
all_crit_met = all_crit_met & ~any_sub0;

%% ============================================================================================== %%
%% HOT-COLD RATIO SELECTION
% see ref. [2] for the criteria
% determine which hc-gain ratios have at least one cal gain scale where criteria are met: 
hc_crit_met = ceil(mean(all_crit_met));
hc_crit_met_inds = find(hc_crit_met);

% convert the gain mults < 1 to being > 1 for evaluating closeness to 1 in a fractional sense:
closeness = VF.hot_mults;
closeness(closeness < 1) = 1 ./ closeness(closeness < 1);
closeness = closeness - 1;
% first entry in the closeness rank vector is the closest to not modifying the hc-gain ratio:
[~, closeness_rank] = sort(closeness); 

% determine which hc gain multiplier index we are selecting:
hc_ind_sel = NaN; % init NaN - invalid
i = 1; % index
keep_looking = 1; % boolean
while i <= VF.n_hot_mults && keep_looking
    test_ind = closeness_rank(i); % try closest first, etc...
    if sum(hc_crit_met_inds == test_ind) > 0
        hc_ind_sel = test_ind; 
        keep_looking = 0; % stop search
    end
    i = i+1;
end

valid_point_found = 1;
if isnan(hc_ind_sel)
    warning('no valid hc indecies found. Setting to middle index for plotting purposes');
    hc_ind_sel = (VF.n_hot_mults+1) / 2;
    valid_point_found = 0;
end
% the value that we enhance our recorded measurement of only the hot load gain by:
hot_mult_sel = VF.hot_mults(hc_ind_sel);
disp(['selected multiplier for hot load gain measurement: ', num2str(hot_mult_sel, 3)]);

%% ============================================================================================== %%
%% CALIBRATION GAIN SELECTION
% In the case that there are multiple calibration gains (i.e. the things that scale hot and cold, as
% well as a, together) that give fully valid parameter space, the goal is to return the one which
% has the lowest scan rate

R_sel = R(:, hc_ind_sel); % scan rate at the selected index
all_crit_met_sel =  all_crit_met(:, hc_ind_sel); % the good cal gain indicies at the selected index
R_choice = R_sel;
R_choice(~all_crit_met_sel) = NaN; 
[R_min, Gcal_ind_sel] = min(R_choice); % minimize scan rate to be conservative 
Gcal_sel_1q_dB = VF.GA_addnoise_scale_test_1q_dB(Gcal_ind_sel);

disp(['min scan rate: ', num2str(R_min), ' (a.u.) at cal gain: ', num2str(Gcal_sel_1q_dB), ...
    ' dB (1q)']);

%% ============================================================================================== %%
%% OUTPUT
% abg outputs:
params_select.abg_out.S_c_1q_qta = abg_test_out.S_c_test_1q_qta(Gcal_ind_sel, hc_ind_sel);
params_select.abg_out.rod_partic = abg_test_out.rod_partic_test(Gcal_ind_sel, hc_ind_sel);
params_select.abg_out.G_s_pow = abg_test_out.G_s_test_pow(Gcal_ind_sel, :, hc_ind_sel);
params_select.abg_out.f_IF_Hz = abg_test_out.f_IF_Hz;
% added noise output: 
params_select.N_Hp_1q_qta = VF.N_Hp_test_1q_qta(Gcal_ind_sel, :, hc_ind_sel); 

% boolean that tells us if a successful point was found in gain space
params_select.valid_point_found = valid_point_found;

%% ============================================================================================== %%
%% PLOT PARAMTERS
hmet_NH = 0.95; % fractional height at which to display criterion met bar
hmet_GS = 0.93;
hmet_Sc = 0.91;
hmet_all = 0.89;

cNH = SD.yellow;
cGS = SD.blue;
calpha = SD.neoncarrot;
cSc = SD.red;
c2 = SD.black;
c_all = SD.black;
c_nsub0 = SD.saddlebrown;

val_mark = '.';
val_sz = 15;
sub0_mark = '+';
sub0_sz = 6;
csub0 = SD.mygray;


if plot_params
    col_fct = @(i,imax,col1,col2) col1*(imax-i)/(imax-1) + col2*(i-1)/(imax-1);
    f_lim = [f_IF_min_int_kHz, f_IF_max_int_kHz]*1E-3;
    ccrit = SD.myred;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HEMT+ added noise
    y_min = -5;
    y_max =  max(max(NH(:,:,hc_ind_sel), [], 'all'), y_min+1E-12);

    fig = figure;
    hold on;
    for j = 1:n_GAs
        lw = 0.5; 
        ls = '--';
        if NH_crit_met(j, hc_ind_sel)
            lw = 2;
            ls = '-';
        end
        col = col_fct(j,n_GAs,cNH,c2);
        plot(f_IF_Hz*1E-6, NH(j,:,hc_ind_sel), 'Color', col, 'LineWidth', lw, 'LineStyle', ls);
    end
    plot(f_lim, VF.crit.NHmin_1q*[1,1], 'Color', ccrit);
    
    ylabel('$N_H$ [1q qta]');
    xlabel('$f_\mathrm{IF}$ [MHz] ');
    xlim(f_lim);
    ylim([y_min, y_max]);
    title(['OP ', num2str(OP_ind), ', $M_{G_{A,h}}=', num2str(hot_mult_sel, 3),'$']);
    hold off;
    if save_figs
        saveas(fig, [fig_folder, f, 'OP ', num2str(OP_ind), ' Added Noise Space'], 'jpeg');
        disp('HEMT+ added noise space figure saved successfully');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % squeezing
    y_min = -0.5;
    y_max = 1.5;
    
    fig = figure;
    hold on;
    for j = 1:n_GAs
        col = col_fct(j,n_GAs,cGS,c2);
        lw = 0.5; 
        ls = '--';
        if GS_crit_met(j,hc_ind_sel)
            lw = 2;
            ls = '-';
        end
        plot(f_IF_Hz*1E-6, GS(j,:,hc_ind_sel), 'Color', col, 'LineWidth', lw, 'LineStyle', ls);
    end
    plot(f_lim, VF.crit.Gs_min_1q_pow*[1,1], 'Color', ccrit);
    plot(f_lim, VF.crit.Gs_max_1q_pow*[1,1], 'Color', ccrit);
    
    ylabel('$G_s$');
    xlabel('$f_\mathrm{IF}$ [MHz] ');
    xlim(f_lim);
    ylim([y_min, y_max]);
    title(['OP ', num2str(OP_ind), ', $M_{G_{A,h}}=', num2str(hot_mult_sel, 3),'$']);
    hold off;
    if save_figs
        saveas(fig, [fig_folder, f, 'OP ', num2str(OP_ind), ' Squeezing Space'], 'jpeg');
        disp('squeezing space figure saved successfully');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cavity spectral density
    x_min = VF.GA_addnoise_scale_min_1q_dB;
    x_max = VF.GA_addnoise_scale_max_1q_dB;
    y_min = min(Sc,[],'all');
    y_max = max(Sc,[],'all');   
    Hmet_NH = y_min+(y_max-y_min)*hmet_NH; % absolute height for criterion-met bar
    Hmet_GS = y_min+(y_max-y_min)*hmet_GS;
    Hmet_Sc = y_min+(y_max-y_min)*hmet_Sc;
    Hmet_all = y_min+(y_max-y_min)*hmet_all;

    fig = figure;
    hold on;
    x_dat = VF.GA_addnoise_scale_test_1q_dB;
    plot(x_dat, Sc(:,hc_ind_sel), 'Color', cSc, 'LineWidth', 0.5, ...
        'LineStyle', '--');
    plot(VF.GA_addnoise_scale_test_1q_dB(Sc_crit_met(:,hc_ind_sel)), ...
        Sc(Sc_crit_met(:,hc_ind_sel),hc_ind_sel), 'Color', cSc, 'LineWidth', 2);
    plot([x_min, x_max], VF.crit.Scmin_1q_qta*[1,1], 'Color', ccrit, 'LineWidth', 1);
    plot([x_min, x_max], VF.crit.Scmax_1q_qta*[1,1], 'Color', ccrit, 'LineWidth', 1);
    
    plot(x_dat(NH_crit_met(:,hc_ind_sel)), ones(size(x_dat(NH_crit_met(:,hc_ind_sel))))*Hmet_NH, ...
        val_mark, 'Color', cNH, 'MarkerSize', val_sz);
    plot(x_dat(GS_crit_met(:,hc_ind_sel)), ones(size(x_dat(GS_crit_met(:,hc_ind_sel))))*Hmet_GS, ...
        val_mark, 'Color', cGS, 'MarkerSize', val_sz);
    plot(x_dat(Sc_crit_met(:,hc_ind_sel)), ones(size(x_dat(Sc_crit_met(:,hc_ind_sel))))*Hmet_Sc, ...
        val_mark, 'Color', cSc, 'MarkerSize', val_sz);
    plot(x_dat(all_crit_met(:,hc_ind_sel)), ones(size(x_dat(all_crit_met(:,hc_ind_sel))))*...
        Hmet_all, val_mark, 'Color', c_all, 'MarkerSize', val_sz);

    plot(x_dat(NH_sub0(:,hc_ind_sel)), ones(size(x_dat(NH_sub0(:,hc_ind_sel))))*Hmet_NH, ...
        sub0_mark, 'Color', csub0, 'MarkerSize', sub0_sz);
    plot(x_dat(GS_sub0(:,hc_ind_sel)), ones(size(x_dat(GS_sub0(:,hc_ind_sel))))*Hmet_GS, ...
        sub0_mark, 'Color', csub0, 'MarkerSize', sub0_sz);
    plot(x_dat(Sc_sub0(:,hc_ind_sel)), ones(size(x_dat(Sc_sub0(:,hc_ind_sel))))*Hmet_Sc, ...
        sub0_mark, 'Color', csub0, 'MarkerSize', sub0_sz);
    plot(x_dat(any_sub0(:,hc_ind_sel)), ones(size(x_dat(any_sub0(:,hc_ind_sel))))*Hmet_all, ...
        sub0_mark, 'Color', csub0, 'MarkerSize', sub0_sz);
    
    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    xlabel('$G_{A,1q}^\mathrm{cal}$ [dB]');
    ylabel('$S_{c}$ [1q qta]');
    title(['OP ', num2str(OP_ind), ', $M_{G_{A,h}}=', num2str(hot_mult_sel, 3),'$']);
    hold off;
    if save_figs
        saveas(fig, [fig_folder, f, 'OP ', num2str(OP_ind), ' Cavity PSD Space'], 'jpeg');
        disp('cavity PSD space figure saved successfully');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % visibility
    y_min = 0;
    if sum(~any_sub0(:, hc_ind_sel)) >= 1
        y_max = 1.15*max(vis(~any_sub0(:, hc_ind_sel),:,hc_ind_sel), [], 'all');
    else
        y_max = 1;
    end
        
    fig = figure;
    hold on;
    for j = 1:n_GAs
        col = col_fct(j,n_GAs,calpha,c2);
        lw = 0.5;
        ls = '--';
        if all_crit_met(j, hc_ind_sel)
            lw = 2;
            ls = '-';
        end
        if ~any_sub0(j, hc_ind_sel)
            plot(f_IF_Hz*1E-6, vis(j,:,hc_ind_sel), 'Color', col, 'LineWidth', lw, 'LineStyle', ls);
        end
    end
    ylabel('$\alpha$');
    xlabel('$f_\mathrm{IF}$ [MHz] ');
    xlim(f_lim);
    ylim([y_min, y_max]);
    title(['OP ', num2str(OP_ind), ', $M_{G_{A,h}}=', num2str(hot_mult_sel, 3),'$']);
    hold off;
    if save_figs
        saveas(fig, [fig_folder, f, 'OP ', num2str(OP_ind), ' Visibility Space'], 'jpeg');
        disp('visibility space figure saved successfully');
    end
end

%% ============================================================================================== %%
%% PLOT VALUE FUNCTION
if plot_valfct
    x_min = VF.GA_addnoise_scale_min_1q_dB;
    x_max = VF.GA_addnoise_scale_max_1q_dB;
    y_min = 0; 
    if sum(~any_sub0,'all') >= 1
        y_max = 1.15 * max(R(~any_sub0));
    else
        y_max = 1;
    end
    Hmet_NH = y_min+(y_max-y_min)*hmet_NH; % absolute height for criterion-met bar
    Hmet_GS = y_min+(y_max-y_min)*hmet_GS;
    Hmet_Sc = y_min+(y_max-y_min)*hmet_Sc;
    Hmet_all = y_min+(y_max-y_min)*hmet_all;
    
    fig = figure;
    hold on;
    x_dat = VF.GA_addnoise_scale_test_1q_dB;
    % only plot where no key paramters go sub-zero
    plot(VF.GA_addnoise_scale_test_1q_dB(~any_sub0(:,hc_ind_sel)), ...
        R(~any_sub0(:,hc_ind_sel), hc_ind_sel), 'LineWidth', 1.5, 'Color', c_nsub0); 
    plot(VF.GA_addnoise_scale_test_1q_dB(all_crit_met(:,hc_ind_sel)), ...
        R(all_crit_met(:,hc_ind_sel), hc_ind_sel), 'LineWidth', 3, 'Color', c_all); 
    plot(Gcal_sel_1q_dB*[1,1], [y_min, y_max], '--', 'LineWidth', 0.75, 'Color', SD.mygray);

    
    plot(x_dat(NH_crit_met(:,hc_ind_sel)), ones(size(x_dat(NH_crit_met(:,hc_ind_sel))))*Hmet_NH, ...
        val_mark, 'Color', cNH, 'MarkerSize', val_sz);
    plot(x_dat(GS_crit_met(:,hc_ind_sel)), ones(size(x_dat(GS_crit_met(:,hc_ind_sel))))*Hmet_GS, ...
        val_mark, 'Color', cGS, 'MarkerSize', val_sz);
    plot(x_dat(Sc_crit_met(:,hc_ind_sel)), ones(size(x_dat(Sc_crit_met(:,hc_ind_sel))))*Hmet_Sc, ...
        val_mark, 'Color', cSc, 'MarkerSize', val_sz);
    plot(x_dat(all_crit_met(:,hc_ind_sel)), ones(size(x_dat(all_crit_met(:,hc_ind_sel))))*...
        Hmet_all, val_mark, 'Color', c_all, 'MarkerSize', val_sz);

    plot(x_dat(NH_sub0(:,hc_ind_sel)), ones(size(x_dat(NH_sub0(:,hc_ind_sel))))*Hmet_NH, ...
        sub0_mark, 'Color', csub0, 'MarkerSize', sub0_sz);
    plot(x_dat(GS_sub0(:,hc_ind_sel)), ones(size(x_dat(GS_sub0(:,hc_ind_sel))))*Hmet_GS, ...
        sub0_mark, 'Color', csub0, 'MarkerSize', sub0_sz);
    plot(x_dat(Sc_sub0(:,hc_ind_sel)), ones(size(x_dat(Sc_sub0(:,hc_ind_sel))))*Hmet_Sc, ...
        sub0_mark, 'Color', csub0, 'MarkerSize', sub0_sz);
    plot(x_dat(any_sub0(:,hc_ind_sel)), ones(size(x_dat(any_sub0(:,hc_ind_sel))))*Hmet_all, ...
        sub0_mark, 'Color', csub0, 'MarkerSize', sub0_sz);
    
    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    xlabel('$G_{A,1q}^\mathrm{cal}$ [dB]');
    ylabel('$R$ [A.U.]');
    title(['ValFunc: OP ', num2str(OP_ind), ', $M_{G_{A,h}}=', num2str(hot_mult_sel, 3),'$']);
    hold off;
    if save_figs
        saveas(fig, [fig_folder, f, 'OP ', num2str(OP_ind), ' Value Fct Params'], 'jpeg');
        disp('value function parameters figure saved successfully');
    end
end

%% ============================================================================================== %%
%% PLOT PARAMETER SPACE VALIDITY
if plot_valid_param_space
    plot_map = NaN(size(all_crit_met));
    % make an ordering: 
    plot_map(~NH_crit_met & ~GS_crit_met & ~Sc_crit_met) = 0;
    plot_map(NH_crit_met & ~GS_crit_met & ~Sc_crit_met) = 1;
    plot_map(~NH_crit_met & GS_crit_met & ~Sc_crit_met) = 2;
    plot_map(~NH_crit_met & ~GS_crit_met & Sc_crit_met) = 3;
    plot_map(NH_crit_met & GS_crit_met & ~Sc_crit_met) = 4;
    plot_map(NH_crit_met & ~GS_crit_met & Sc_crit_met) = 5;
    plot_map(~NH_crit_met & GS_crit_met & Sc_crit_met) = 6;
    plot_map(NH_crit_met & GS_crit_met & Sc_crit_met) = 7;
    c_none = SD.white;
    cNHGS = SD.green;
    cNHSc = SD.orange;
    cGSSc = SD.purple;
    c_vec = [c_none; cNH; cGS; cSc; cNHGS; cNHSc; cGSSc; c_all];
    
    
    hm_step = VF.hot_mults(2) - VF.hot_mults(1);
    GA_step = VF.GA_addnoise_scale_test_1q_dB(2) - VF.GA_addnoise_scale_test_1q_dB(1);
    plot_hot_mults = [VF.hot_mults, VF.hot_mults(end) + hm_step];
    plot_hot_mults = plot_hot_mults - hm_step/2;
    plot_GA = [VF.GA_addnoise_scale_test_1q_dB, VF.GA_addnoise_scale_test_1q_dB(end) + GA_step];
    plot_GA = plot_GA - GA_step/2;
    plot_map = [plot_map, NaN(size(plot_map,1),1)];
    plot_map = [plot_map; NaN(1, size(plot_map,2))];
    
    % these next two lines are important to getting the color scaling right. They do not display in
    % the plot because of how pcolor works
    plot_map(1, end) = 0; 
    plot_map(end, 1) = 7;
    
    x_min = min(VF.hot_mults) - hm_step/2;
    y_min = min(VF.GA_addnoise_scale_test_1q_dB) - GA_step/2;
    x_max = max(VF.hot_mults) + hm_step/2;
    y_max = max(VF.GA_addnoise_scale_test_1q_dB) + GA_step/2;
    
    fig = figure;
    hold on;
    pcolor(plot_hot_mults, plot_GA, plot_map);
    colormap(c_vec);
    
    % line at hc-multiplier = 1, the preferred value:
    plot([1,1], [y_min,y_max], '-', 'LineWidth', 1.5);
    % b&w crosshairs at the selected value:
    plot(VF.hot_mults(hc_ind_sel)*[1,1], [y_min,y_max], '-', 'LineWidth', 0.75, 'Color', SD.white);
    plot(VF.hot_mults(hc_ind_sel)*[1,1], [y_min,y_max], '--', 'LineWidth', 0.75, 'Color', SD.black);
    plot([x_min, x_max], Gcal_sel_1q_dB*[1,1], '-', 'LineWidth', 0.75, 'Color', SD.white);
    plot([x_min, x_max], Gcal_sel_1q_dB*[1,1], '--', 'LineWidth', 0.75, 'Color', SD.black);

    
    % fill frame borders back in:
    plot([x_min,x_max], [y_min,y_min]);
    plot([x_min,x_max], [y_max,y_max]);
    plot([x_min,x_min], [y_min,y_max]);
    plot([x_max,x_max], [y_min,y_max]);
    
    xlim([x_min,x_max]);
    ylim([y_min,y_max]);
    xlabel('$M_{G_{A,h}}$');
    ylabel('$G^\mathrm{cal}_{1q}$');
    title(['OP ', num2str(OP_ind)]);
    hold off;
    
    if save_figs
        saveas(fig, [fig_folder, f, 'OP ', num2str(OP_ind), ' Valid Param Space'], 'jpeg');
        disp('valid parameter space figure saved successfully');
    end

end

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################