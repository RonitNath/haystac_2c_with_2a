%% DESCRIPTION
% processing code for HAYSTAC phase 2 data
%
%% HISTORY
% - v1 created by Dan Palken on 8/27/19. Based directly off of HAYSTACp2_proto_process_v11 code
% - v2 created by DP on 8/28/19. Goes with wkspc_construction_v1
% - v3 created by DP on 8/30/19
% - v4 created by DP on 8/30/19
% - v5 created by DP on 9/3/19. Allows for different frequency vectors for the different ENA
% measurements
% - v6 created by DP on 9/3/19
% - v7 created by DP on 9/4/19
% - v8 created by DP on 9/5/19. Plots all raw spectra together
% - v9 created by DP on 9/5/19. Adds a new step of multiplying all the sraw spectra to have the same
% mean before IF cuts are made
% - v10 created by DP on 9/8/19
% - v11 created by DP on 9/9/19. Incorporates BPM exclusion based on ref. [22]
% - v12 created by DP on 9/10/19. Reintroduces NaN's into the combined spectrum. Smarter NaN summing
% - v13 created by DP on 9/13/19
% - v14 created by DP on 9/16/19
% - v15 created by DP on 9/16/19. Better IF mean excess plotting
% - v16 created by DP on 9/16/19
% - v17 created by DP on 9/17/19. Plots prior updates
% - v18 created by DP on 9/19/19
% - v19 created by DP on 9/20/19. Incorporates BPM exclusion based on ref. [24]
% - v20 created by DP on 9/24/19. Plots cavity tuning steps
% - v21 created by DP on 9/24/19. Plots cavity tuning steps
% - v22 created by DP on 11/10/19. Works with wkspc_construction_v9 outputs
% - v23 created by DP on 1/22/20. Includes code to make manual rf
% - v24 created by DP on 1/23/20. Moves manual rf cuts to after all of the SG filtering
% - v25 created by DP on 1/25/20. Combines the processing of many datasets
% - v26 created by DP on 2/16/20
% - v27 created by DP on 2/16/20. Breaks IF processing into batches
% - v28 created by DP on 2/25/20
% - v29 created by DP on 3/1/20. 3 batches, spectrum cuts
% - v30 created by DP on 3/6/20. Different IF windows for different batches
% - v30a created by DP on 3/6/20. Alternate attempt at making different IF bands from v30. Lower Kg
% - v31 created by DP on 3/18/20. Aligns a bit more with Kelly's procedure
% - v31 created by DP on 3/19/20. Fixed AMP gains to proper power units
% - v32 created by DP, date uncertain
% - v33 created by DP on 3/26/20. Outputs exlcusion data for rescan processing code
% - v34 created by DP on 4/11/20
% - v35 created by DP on 4/14/20
%

%% REFERENCES
% [1]: DP's lab NB: "Data Acquisition + Analysis > Mock-Axion Experiment > Processing the
% Spectra/Disambiguation," 4 Mar 2018
% [2]: B. Brubaker et al, Phys Rev D. 96, 123008 (2017)
% [3]: DP's lab NB: "Data Acquisition + Analysis > Mock-Axion Experiment > Processing the
% Spectra/Processing the Spectra, part 2"
% [4]: MM's lab NB: "measurements > scan rate enhancment > cold cavity spectro," top of page
% [5]: DP's lab NB: "Data Acquisition + Analysis > Mock-Axion Experiment > Measuring the SNR >
% Theoretical Alternative," 2 May 2018
% [6]: DP's lab NB: "Data Acquisition + Analysis > Mock-Axion Experiment > Measuring the SNR >
% Refined Model and Method," 2 May 2018
% [7]: DP's lab NB: "Data Acquisition + Analysis > Mock-Axion Experiment > Noise Floor Data,"
% 4-7 May 2018
% [8]: MATLAB code: noise_profile_analysis_v1
% [9]: DP's lab NB: "Data Acquisition + Analysis > Mock-Axion Experiment > Close Reading of Ben's
% PRD > PRD Reading, part 3," 15 May 2018-19 Jun 2018
% [10]: DP's lab NB: "Data Acquisition + Analysis > Mock-Axion Experiment > Processing the
% Spectra/Processing the Spectra, part 6," 28 Jun 2018
% [11]: DP's lab NB: "Data Acquisition + Analysis > Mock-Axion Experiment > Close Reading of Ben's
% PRD > PRD Reading, part 4," 29 Jun 2018
% [12]: MATLAB 2018b documentation: Chi-Square Goodness-of-Fit Test
% [13]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > HAYSTACp2 Processing Running
% Jounral > Axion Lineshape"
% [14]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > HAYSTACp2 Processing Running
% Jounral > Misalignment"
% [15]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > HAYSTACp2 Processing Running
% Jounral"
% [16]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > HAYSTACp2 Processing Running
% Jounral > Filling the Grand Spectrum"
% [17]: DP's lab NB: "Data Acquisition + Analysis > Mock-Axion Experiment > Processing the
% Spectra/Processing the Spectra, part 3"
% [18]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Squeezing/Hot Rod
% Calibration > Measurmenet Instructions for Kelly Backes"
% [19]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Live Data," 5 Sep 2019
% [20]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Statistical Testing for
% Gaussianity"
% [21]: "A Maximally informative axion haloscope analysis," manuscript, Dan Palken
% [22]: MATLAB code: Bayesian_Analyze_v10
% [23]: DP's lab NB, 'Publications/A maximally informative axion haloscope analysis/Intelligent
% Smoothing/Averaging,' 9/12/18
% [24]: MATLAB code: Bayesian_Analyze_v16
% [25]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Spectrum Cuts"
%

%% ============================================================================================== %%
%% INITIALIZATION
clear;
close all;
addpath('MATLAB Universal')
delete(findall(0,'type','figure','tag','TMWWaitbar'));
dt_start = datetime; % record start time
Settings_and_Defaults_v4; % apply custom settings
ws = warning('on', 'all');
set(groot, 'defaultFigureColorMap', flipud(gray(10000))); % smooth out the color scheme
% set(0, 'DefaultFigureWindowStyle', 'docked'); % docks the plots
set(0, 'DefaultFigureWindowStyle', 'normal'); % undocks the plots


%% ============================================================================================== %%
%% CONSTANTS
% fundamnetal constants:
kB_JperK = 1.38064852E-23; % Boltzmann constant;
e_C = 1.6021766208E-19; % elementary charge
h_Planck_Js = 6.626070040E-34; % Planck's constant
hbar_Js = h_Planck_Js / (2 * pi);% reduced Planck's constant
c_mpers = 2.99792458E8; % speed of light
mu0_Hperm = 4*pi*1E-7; % permeability of free space
epsilon0_Fperm = 1 / (c_mpers^2 * mu0_Hperm); % permittivity of free space
alpha_fs = 1 / (4 * pi * epsilon0_Fperm) * e_C^2 / (hbar_Js * c_mpers); % fine structure constant
vc_mpers = 220E3; % mode of the axion Maxwellian velocity distribution in the galactic rest frame
Ev2_m2pers2 = (3/2) * vc_mpers^2; % see refs. [2], [13]
EBeta2 = Ev2_m2pers2 / c_mpers^2; % see refs. [2], [13]
vs_mpers = 232E3; % velocity of sun about galactic center; see refs. [2], [13]
r = vs_mpers / sqrt(Ev2_m2pers2); % will be sqrt(2/3) when vs = vc; see refs. [2], [13]

tab = '    '; % for dispaying messages
f = filesep; % shortcut for platform-specific file separator

% putative axion constants:
g_gamma_KSVZ = 0.97;
rho_a_GeVpercm3 = 0.45; % local dark matter energy density
rho_a_Jperm3 = rho_a_GeVpercm3 * e_C*1E9*1E6; % in SI units
Lambda_MeV = 75.6; % encodes dependance on hadronic physics. Ben endorses this value
Lambda_J = Lambda_MeV * e_C*1E6; % SI units

%% ============================================================================================== %%
%% USER INPUT - PROCESSING + PLOTTING
% general user inputs or those for processing the data

ver = 'v35'; % should match script name. Not for the data to be loaded

wkspc_construction_ver = 'v16';
% names of the data runs
data_runs = {}; % init empty
run_batch = [];
% add runs one-by-one for ease of commenting out
data_runs = [data_runs, '20221111']; run_batch = [run_batch, 1];
% data_runs = [data_runs, '20190903']; run_batch = [run_batch, 1];
% data_runs = [data_runs, '20190907']; run_batch = [run_batch, 1];
% data_runs = [data_runs, '20190911']; run_batch = [run_batch, 1];
% data_runs = [data_runs, '20190915']; run_batch = [run_batch, 1];
% data_runs = [data_runs, '20190917']; run_batch = [run_batch, 1];
% data_runs = [data_runs, '20190919']; run_batch = [run_batch, 1]; 
% data_runs = [data_runs, '20190925']; run_batch = [run_batch, 2]; % begin batch 2
% data_runs = [data_runs, '20190930']; run_batch = [run_batch, 2];
% data_runs = [data_runs, '20191002']; run_batch = [run_batch, 2];
% data_runs = [data_runs, '20191003']; run_batch = [run_batch, 3]; % begin batch 3
% data_runs = [data_runs, '20191008']; run_batch = [run_batch, 3];
% data_runs = [data_runs, '20191011']; run_batch = [run_batch, 3];
% data_runs = [data_runs, '20191014']; run_batch = [run_batch, 3];
% data_runs = [data_runs, '20191019']; run_batch = [run_batch, 3];
% data_runs = [data_runs, '20191023']; run_batch = [run_batch, 3];
% data_runs = [data_runs, '20191028']; run_batch = [run_batch, 3];
% data_runs = [data_runs, '20191102']; run_batch = [run_batch, 3];
% data_runs = [data_runs, '20191107']; run_batch = [run_batch, 3];
% data_runs = [data_runs, '20191111']; run_batch = [run_batch, 3];
% data_runs = [data_runs, '20191125']; run_batch = [run_batch, 3];
% data_runs = [data_runs, '20191130']; run_batch = [run_batch, 3];
% data_runs = [data_runs, '20191205']; run_batch = [run_batch, 3];
% data_runs = [data_runs, '20191209']; run_batch = [run_batch, 3];
% data_runs = [data_runs, '20191214']; run_batch = [run_batch, 3];

n_batches = max(run_batch);

% manual rf cuts
% N_man manual cuts are specified in the following form: an N_man x 2 matrix specifying the start
% and stop rf frequencies to cut at, alongside an N_man x 2 vector of OP indices indicating to which
% spectra (start and end, inlcusive) the cuts apply (this allows exemptions) - it is fine to list a
% range broader than the last spectra which technically encompass data with a cut
hi = 1000000;
f_man_rf_cut_GHz = [4.11062 - 0.0001, 4.11062 + 0.0001; ...
    4.16526 - 0.0001, 4.16526 + 0.0001; ...
    4.158275 - 0.0001, 4.158275 + 0.0001; ...
    4.11471 - 0.0001, 4.11471 + 0.0001; ...
    4.16911 - 0.0001, 4.16911 + 0.0001;];
man_cut_indices = [1,hi; ...
    1,hi; ...
    1,hi; ...
    1,hi; ...
    1,hi; ...
    1,hi];

% booleans:
print_timing = 01;
rand_indiv_OP = 0; % select a random operating point to plot individual spectrum data from
% shortcutting and bypassing:
shortcut_IF_cuts = 01; % if 1, uses previously saved IF cuts. Auto-disables if no file found
shortcut_SG_filtering = 01; % uses previously saved SG filters. Auto-disables if no file found
save_IF_cuts_to_wkspc = 01; % only applies if shortcut_IF_cuts is 0
save_SG_filters_to_wkspc = 01; % only applies if shortcut_SG_filtering is 0
save_scan1excl_to_wkspc = 01; 
bypass_cav_spec = 1; % bypasses all cavity spectroscopy not at the selected operating point
% plotting:
plot_input_params = 01;
plot_raw_spectra_IF = 00;
plot_spec_cuts = 00;
plot_raw_spectrum_RF = 00;
plot_various_mean_spectra_IF = 00;
plot_normalized_excess_RF = 00;
plot_hist_mean_excess = 00;
plot_n_cuts = 00;
plot_proc_spectrum_RF = 00;
plot_waterfall_spectra_RF = 00;
plot_hist_all_proc_spec = 00;
plot_cav_spec = 00;
plot_noise_profile = 00;
plot_resc_SNR = 00;
plot_resc_spectrum_RF = 00;
plot_comb_spectrum = 00;
plot_and_hist_norm_comb_spectrum = 010
plot_comb_autocor = 01;
plot_axion_PDF = 01;
plot_Lq = 01;
plot_and_hist_norm_grand_spectrum = 01;
plot_grand_autocor = 01;
plot_and_hist_renorm_grand_spectrum = 01;
plot_excl = 01;
save_figs = 01;
% for histogramming some pre-comibined spectra, plot a subset of the naive and precise expectations:
show_precomb_expected_Gaussians = 01;
% for analysis plotting:
show_LU_inset = 01;
show_speedup_inset = 01;
mark_rescan_freqs = 01;

indiv_OP_ind = 702; % iff random selection is off, use this OP for individual spectrum data

% multiply the raw spectra by some amount at the time of loading to avoid potential numerical
% precision errors later on:
raw_initial_multiplier = 1E12; % set to 1 to leave unchanged. 1E12 seems about right for HAYSTACp2

% processing band: go about a mumber of FWHM's in either direction from the center of the cavity
max_processing_freq_IF_MHz = 1.34505;
max_analysis_freq_IF_MHz = [0.74505, 1.34505, 1.34505];
min_analysis_freq_IF_MHz = 0.04505; % same for all

% filtering:
freq_domain_SG = 0; % apply the SG filter as a fit in frequency domain. 2 is native functionality
% parameters for the filter applied to the IF-averaged raw spectra:
SG_raw_filt_deg = 10;
SG_raw_filt_W_kHz = 50; % "half-width" of window (full actually 2W + 1)
% parameters for the main filter applied later on:
SG_norm_filt_deg = 4;
SG_norm_filt_W_kHz = SG_raw_filt_W_kHz;

% flagging for cuts:
n_adj_flag_IF = 3;
n_adj_flag_RF = 3;

RF_cutoff = 6; % units of proc. spec. sigmas for waterfall plotting

% IF cuts (percentage of cut bins may not match expectation if not in Gaussian limit):
% FLAG - at lower values, e,g, 4.5 in v4, this causes a crash
IF_discard_thresh_sigma = 4.5;
% RF cuts (percentage of cut bins may not match expectation if not in Gaussian limit):
RF_discard_thresh_sigma = 6;

% axion PDF
ax_int_show_det_min_kHz = -1.5; % 0 = start at the rest mass
ax_int_kHz = 9.0; % integrate axion this many linewidths out - multiple of raw spectum bin spacing
grand_density = 1; % have integer no. of grand bins per combined bin (opposite of ref. [2])

corr_lag_extent_MHz = 1; % optionally look at correlations this far out

% Kelly's simulations give 0.91 for eta instead of 0.90. For rescans she gets 0.90
eta_p = 0.91 * 0.93; % see p. 21 of ref. [2]. % FLAG - update to 0.90 for rescans

% save-to:
n_runs = length(data_runs);
if n_runs == 1
    runs_name = data_runs{1};
else
    runs_name = 'multirun';
end
save_to_folder = ['save_to_', runs_name, '_proc_', ver];

% for shortcutting:
save_to_IF_cuts = 'IF_cuts.mat';
save_to_SG_filters = 'SG_filters.mat';
% for rescan processing: 
save_to_scan1excl = 'scan1excl.mat';

%% ============================================================================================== %%
%% USER INPUT - ANALYSIS
% user inputs for Bayesian Power Measured (BPM) and Frequentist Threshold (FT) analysis of the data,
% ref. [21]
% boolean(s):
norm_cand_mks_szs_to_priors = 0;

% consider many couplings (normalized to KSVZ)
g_target_min = 0.5;
g_target_max = 2.0;
g_target_step = 0.01;
g_target_show_top = g_target_max; % only show up to an interesting coupling (<= g_target_max)

u_flr = 0.05;
u_ceil = 1 / u_flr;
desired_prior_red = 0.1; % e.g. 0.1 ~ "90% exlcusion"

% averaging the data to have many fewer frequencies to plot (~200,000 is too many to see)
n_red_freqs = 100;

% set lower and upper cutoff frequencies, here as percentages, for aggregation:
% FLAG - very narrow range
low_cutoff_freq_pct = 1; % expressed as a number 0-100
hi_cutoff_freq_pct = 99; % expressed as a number 0-100

% identify a specified number of "candidates" within the range, above
n_cands = 10;

% frequentism:
mu_target = 5.1; % look for 5.1 sigma axions in line with phase 1, run 1 (somewhat arbitrary)
n_thresh_scans = 2; % note this explicitely assumed to be 2 at times

% rescans: 
rescan_cluster_within_Hz = 250E3; 

%% ============================================================================================== %%
%% DATA RUN INFO
% directory and file names to load raw spectra and accompanying info from:
OP_dat_name = 'dat_HAYSTACp2_OP';

% now load in some infomration needed for allocation:
n_OPs = 0; % init 0
for n = 1:n_runs
    OP_dir = ['HAYSTACp2_raw_saved_', data_runs{n}, '_', wkspc_construction_ver];
    c = 0; % init 0 for each data run
    while exist([OP_dir, f, OP_dat_name, num2str(c+1), '.mat'], 'file')
        c = c + 1;
    end
    n_OPs = n_OPs + c;
end
disp([num2str(n_OPs), ' opearting point workspaces detected']);
man_cut_indices(man_cut_indices > n_OPs) = n_OPs; % update manual cut indecies in light of n_OPs

% choose a directory to use as the temporary one:
% this should be chosen such that it is the one with the most IF bins
OP_dir_temp = ['HAYSTACp2_raw_saved_', data_runs{length(data_runs)}, '_', wkspc_construction_ver];

temp_first_OP = load([OP_dir_temp, f, OP_dat_name , num2str(1)]).OP_data;
n_ENA_fs_tx1 = length(temp_first_OP.ENA_f_tx1_GHz);
n_ENA_fs_tx2 = length(temp_first_OP.ENA_f_tx2_GHz);
n_ENA_fs_rf = length(temp_first_OP.ENA_f_rf_GHz);
n_ENA_fs_G_AMP = length(temp_first_OP.ENA_f_G_AMP_GHz);
n_fs4add_noise = length(temp_first_OP.f4NH_IF_Hz);
n_fs4abg = length(temp_first_OP.abg_out.f_IF_Hz);
n_IF_fs = length(temp_first_OP.f_IF_MHz);
raw_res_Hz = round(1E6 * (temp_first_OP.f_IF_MHz(3) - temp_first_OP.f_IF_MHz(2)),3);
clear(varname(temp_first_OP));

%% ============================================================================================== %%
%% ALLOCATE INTAKE
% allocate only those inputs which may or may not vary by operating point:
intake.batch = NaN(n_OPs, 1);
intake.sub_OP_num = NaN(n_OPs, 1);
intake.tone_height = NaN(n_OPs, 1);
intake.tone_height_std = NaN(n_OPs, 1);
intake.BW_AMP_MHz = NaN(n_OPs, 1);
intake.f_cav_tx1_GHz = NaN(n_OPs, 1);
intake.f_cav_tx2_GHz = NaN(n_OPs, 1);
intake.f_cav_r1_GHz = NaN(n_OPs, 1);
intake.f_cav_r2_GHz = NaN(n_OPs, 1);
intake.G_AMP_1q_pow = NaN(n_OPs, 1);
intake.direct_squeezing_pow = NaN(n_OPs, 1);
intake.cav_kappa_ext_tx1_kHz = NaN(n_OPs, 1);
intake.cav_kappa_loss_tx1_kHz = NaN(n_OPs, 1);
intake.cav_kappa_ext_tx2_kHz = NaN(n_OPs, 1);
intake.cav_kappa_loss_tx2_kHz = NaN(n_OPs, 1);
intake.C_010 = NaN(n_OPs, 1);
intake.T_fridge_K = NaN(n_OPs, 1);
intake.S_c_1q_qta = NaN(n_OPs, 1);
intake.rod_partic = NaN(n_OPs, 1);
intake.G_s_pow = NaN(n_OPs, n_fs4abg);
intake.N_Hp_1q_qta = NaN(n_OPs, n_fs4add_noise);
intake.ENA_f_tx1_GHz = NaN(n_OPs, n_ENA_fs_tx1);
intake.ENA_f_tx2_GHz = NaN(n_OPs, n_ENA_fs_tx2);
intake.ENA_f_rf_GHz = NaN(n_OPs, n_ENA_fs_rf);
intake.ENA_f_G_AMP_GHz = NaN(n_OPs, n_ENA_fs_G_AMP);
intake.ENA_Ssw_tx1 = NaN(n_OPs, n_ENA_fs_tx1);
intake.ENA_Ssw_tx2 = NaN(n_OPs, n_ENA_fs_tx2);
intake.ENA_Sss_rf = NaN(n_OPs, n_ENA_fs_rf);
intake.raw_spec_IF = NaN(n_OPs, n_IF_fs);

%% ============================================================================================== %%
%% LOAD OP DATA
% load all operating point data into one big structure of matrices and vectors
disp('loading haloscope operational and spectral data...');
wb = waitbar(0, 'loading OP data...');

abs_i = 1; % continues counting from run to run
for n = 1:n_runs
    OP_dir = ['HAYSTACp2_raw_saved_', data_runs{n}, '_', wkspc_construction_ver];
    i = 1; % resets to 1 at each data run
    while exist([OP_dir, f, OP_dat_name, num2str(i), '.mat'], 'file') % same condition as above
        load([OP_dir, f, OP_dat_name, num2str(i)]);
         % inputs which are guaranteed to be the same for all operating points - make sure to choose
         % operating point with most IF bins
        if abs_i == n_OPs
            intake.f_IF_MHz = OP_data.f_IF_MHz;
            intake.B0_T = OP_data.B0_T;
            intake.V_m3 = OP_data.V_m3;
            intake.sub_per_raw = OP_data.sub_per_raw;
            intake.SQ_cav_eta = OP_data.SQ_cav_eta;
            intake.cav_AMP_eta = OP_data.cav_AMP_eta;
            intake.W_width_MHz = OP_data.W_width_MHz;
            intake.W_center_MHz = OP_data.W_center_MHz;
            intake.T_rod_K = OP_data.T_rod_K;
            intake.f_abg_Hz = OP_data.abg_out.f_IF_Hz;
            intake.f_NH_Hz = OP_data.f4NH_IF_Hz;
        end
        intake.date{abs_i} = data_runs{n};
        intake.sub_OP_num(abs_i) = i;
        intake.batch(abs_i) = run_batch(n);
        intake.ENA_f_tx1_GHz(abs_i, :) = OP_data.ENA_f_tx1_GHz;
        intake.ENA_f_tx2_GHz(abs_i, :) = OP_data.ENA_f_tx2_GHz;
        intake.ENA_f_rf_GHz(abs_i, :) = OP_data.ENA_f_rf_GHz;
        intake.ENA_f_G_AMP_GHz(abs_i, :) = OP_data.ENA_f_G_AMP_GHz;
        intake.tone_height(abs_i) = OP_data.tone_height;
        intake.tone_height_std(abs_i) = OP_data.tone_height_std;
        intake.BW_AMP_MHz(abs_i) = OP_data.BW_AMP_MHz;
        intake.f_cav_tx1_GHz(abs_i) = OP_data.f_cav_t1_GHz;
        intake.f_cav_tx2_GHz(abs_i) = OP_data.f_cav_t2_GHz;
        intake.f_cav_r1_GHz(abs_i) = OP_data.f_cav_r_GHz;
        intake.f_cav_r2_GHz(abs_i) = OP_data.f_cav_r2_GHz;
        intake.G_AMP_1q_pow(abs_i) = OP_data.G_AMP_1q_pow;
        intake.direct_squeezing_pow(abs_i) = OP_data.direct_squeezing;
        intake.S_c_1q_qta(abs_i) = OP_data.abg_out.S_c_1q_qta;
        intake.rod_partic(abs_i) = OP_data.abg_out.rod_partic;
        intake.G_s_pow(abs_i, :) = OP_data.abg_out.G_s_pow;
        % FLAG - make sure code knows it always uses 1q:
        intake.N_Hp_1q_qta(abs_i, :) = OP_data.N_Hp_1q_qta;
        intake.T_fridge_K(abs_i) = OP_data.T_fridge_K;
        intake.cav_kappa_ext_tx1_kHz(abs_i) = OP_data.cav_kappa_ext_t1_kHz;
        intake.cav_kappa_loss_tx1_kHz(abs_i) = OP_data.cav_kappa_loss_t1_kHz;
        intake.cav_kappa_ext_tx2_kHz(abs_i) = OP_data.cav_kappa_ext_t2_kHz;
        intake.cav_kappa_loss_tx2_kHz(abs_i) = OP_data.cav_kappa_loss_t2_kHz;
        intake.C_010(abs_i) = OP_data.C_010;
        intake.ENA_Ssw_tx1(abs_i, :) = OP_data.ENA_Ssw_t1;
        intake.ENA_Ssw_tx2(abs_i, :) = OP_data.ENA_Ssw_t2;
        intake.ENA_Sss_rf(abs_i, :) = OP_data.ENA_Sss;
        
        intake.raw_spec_IF(abs_i, :) = ...
            [OP_data.raw_spec_IF; NaN(n_IF_fs-length(OP_data.raw_spec_IF),1)];
        
        clear OP_data; % clear to make sure a new one is loaded each time
        
        i = i + 1;
        abs_i = abs_i + 1;
        waitbar(abs_i/n_OPs, wb);
    end
end
close(wb);
disp([tab, 'OP data loaded for ', num2str(n_runs), ' runs']);

%% ============================================================================================== %%
%% DERIVED QUANTITIES - PROCESSING
n_OPs_per_batch = NaN(1, n_batches);
for i = 1:n_batches
    n_OPs_per_batch(i) = sum(intake.batch == i);
end

min_f_IF_orig_MHZ = min(intake.f_IF_MHz);
max_f_IF_orig_MHZ = max(intake.f_IF_MHz);
% rounding gets rid of numerical error
grand_spacing_GHz = raw_res_Hz*1E-9 / grand_density;
Kg = round(ax_int_kHz*1E3 / raw_res_Hz) + 1; % rounding to rid numerical error (already an integer)
% see ref. [15], 5/17/19 for the + 1 (basically Kg is how many CS bins contribute to each GS bin)
if mod(grand_density, 1) ~= 0
    error('non-integral grand spectrum density: ML weighting not designed to accomodate');
end
if mod(ax_int_kHz*1E3 / raw_res_Hz, 1) ~= 0
    error('axion integration window is non-intreger nyumber of raw spectrum bins');
end

f_LO_GHz = intake.f_cav_tx1_GHz; % first measured cavity frequecy is the LO frequency
% compute best estimates cavity frequency & couplings as average of before/after measurements:
f_cav_avg_GHz = (intake.f_cav_tx1_GHz + intake.f_cav_tx2_GHz) / 2;
cav_kappa_ext_avg_kHz = (intake.cav_kappa_ext_tx1_kHz + intake.cav_kappa_ext_tx2_kHz) / 2;
cav_kappa_loss_avg_kHz = (intake.cav_kappa_loss_tx1_kHz + intake.cav_kappa_loss_tx2_kHz) / 2;
cav_kappa_tot_avg_kHz = cav_kappa_ext_avg_kHz + cav_kappa_loss_avg_kHz;
cav_kappa_diff_kHz = cav_kappa_ext_avg_kHz - cav_kappa_loss_avg_kHz; % positive ~ overcoupled

raw_IF_mult = intake.raw_spec_IF * raw_initial_multiplier;
% present math ensures that the number of IF analysis bins will be integral

% bins that get NaN's in the RF spectra because they are below the (noisy) low-end IF cutoff:
n_noisy_center_bins = round(2*1E6 * min_f_IF_orig_MHZ / raw_res_Hz - 1); % will be odd
n_raw_bins_IF = length(intake.f_IF_MHz);
n_raw_bins_RF = 2 * n_raw_bins_IF + n_noisy_center_bins;

SG_raw_filt_W_bins = (SG_raw_filt_W_kHz / raw_res_Hz)*1E3;
SG_norm_filt_W_bins = (SG_norm_filt_W_kHz / raw_res_Hz)*1E3;

% choose a random or pre-specified spectrum index for the plotting of all individual spectra
if rand_indiv_OP % choose an OP randomly
    disp('random individual OP selection ON');
    indiv_OP_ind = randi(n_OPs);
else % stick with the pre-selected value
    disp('random individual OP selection OFF'); %#ok<*UNRCH>
end
disp([tab 'using indvidual OP no. ', num2str(indiv_OP_ind)]);


% axion PDF:
% define a "typical" axion frequency that will be used to calculate all the Lqs:
ax_f_typ_GHz = mean(f_LO_GHz);
ax_frac_mean = EBeta2 * (1 + r^2)/2; % see ref [13]
ax_frac_lw = EBeta2 * sqrt(r^2/3 + 1/6); % see ref. [13] (general case, linewidth := sigma)
ax_mean_GHz = ax_f_typ_GHz + ax_frac_mean * ax_f_typ_GHz;
ax_lw_kHz = ax_frac_lw * ax_f_typ_GHz*1E6;
ax_int_min_GHz = ax_f_typ_GHz + ax_int_show_det_min_kHz*1E-6;
ax_int_max_GHz = ax_f_typ_GHz + ax_int_kHz*1E-6;

% processing & analysis band manipulations
max_processing_freq_p_buff_IF_MHz = max_processing_freq_IF_MHz + ...
    max(SG_norm_filt_W_kHz, SG_raw_filt_W_kHz)*1E-3;
n_processing_bins_RF = round(2 * max_processing_freq_IF_MHz*1E6 / raw_res_Hz);
n_processing_bins_IF = (n_processing_bins_RF - n_noisy_center_bins) / 2;
n_processing_p_buff_bins_IF = n_processing_bins_IF + ...
    max(SG_norm_filt_W_bins, SG_raw_filt_W_bins);
n_processing_p_buff_bins_RF = 2 * n_processing_p_buff_bins_IF + n_noisy_center_bins;
folded_bin_range = (n_processing_bins_RF - n_processing_bins_IF + 1):n_processing_bins_RF;
n_analysis_bins_RF = round(2 * max_analysis_freq_IF_MHz*1E6 / raw_res_Hz);
n_analysis_bins_IF = (n_analysis_bins_RF - n_noisy_center_bins) / 2;

% manual rf cuts:
N_man = size(f_man_rf_cut_GHz, 1);

% HAYSTAC power-conversion constant from ref. [2], Eq. (6):
U0_J = (g_gamma_KSVZ^2 * alpha_fs^2 * hbar_Js^3 * c_mpers^3 * rho_a_Jperm3 * 2 * pi * ...
    intake.cav_AMP_eta * intake.B0_T^2 * intake.V_m3) / (pi^2 * Lambda_J^4 * mu0_Hperm);

% cevity drift
cav_drift_tx_kHz = 1E6 * (intake.f_cav_tx2_GHz - intake.f_cav_tx1_GHz);
cav_drift_rfl_kHz = 1E6 * (intake.f_cav_r2_GHz - intake.f_cav_r1_GHz); 

% create directory if necessary:
if ~exist(save_to_folder, 'dir')
    mkdir(save_to_folder);
end

%% ============================================================================================== %%
%% PLOT INTAKE PARAMS
% plot the intake parameters vs. frequency
if plot_input_params
    h_mult = 1.25;
    w_mult = 1.075;
    b_sub = 0.04;
    pos_adj = [1,1-b_sub,w_mult,h_mult];
    fig = figure;
    
    min_x = 0.5;
    max_x = n_OPs + 0.5;
    
    % cavity parameters:
    subp1 = subplot(4, 2, 1); % cavity frequency
    set(subp1, 'Position', get(subp1, 'Position') .* pos_adj);
    % title('Cavity Intake');
    hold on;
    plot(1:n_OPs, intake.f_cav_tx1_GHz, '.', 'Color', SD.myblue, 'MarkerSize', 25);
    plot(1:n_OPs, intake.f_cav_tx2_GHz, '.', 'Color', SD.myred, 'MarkerSize', 12);
    ylabel('$f_\mathrm{cav}$ (GHz)');
    mylegend({'pre-', 'post-'}, 'Location', 'SouthEast');
    xlim([min_x, max_x]);
    ylim([min([intake.f_cav_tx1_GHz; intake.f_cav_tx2_GHz]), ...
        max([intake.f_cav_tx1_GHz; intake.f_cav_tx2_GHz])]);
    set(gca,'XTick',[]);
    hold off;
    
    subp2 = subplot(4, 2, 3); % measurement port coupling
    set(subp2, 'Position', get(subp2, 'Position') .* pos_adj);
    hold on;
    plot(1:n_OPs, intake.cav_kappa_ext_tx1_kHz, '.', 'Color', SD.myblue, 'MarkerSize', 25);
    plot(1:n_OPs, intake.cav_kappa_ext_tx2_kHz, '.', 'Color', SD.myred, 'MarkerSize', 12);
    ylabel('$\kappa_m$ (kHz)');
    xlim([min_x, max_x]);
    set(gca,'XTick',[]);
    hold off;
    
    subp3 = subplot(4, 2, 5); % loss port coupling
    set(subp3, 'Position', get(subp3, 'Position') .* pos_adj);
    hold on;
    plot(1:n_OPs, intake.cav_kappa_loss_tx1_kHz, '.', 'Color', SD.myblue, 'MarkerSize', 25);
    plot(1:n_OPs, intake.cav_kappa_loss_tx2_kHz, '.', 'Color', SD.myred, 'MarkerSize', 12);
    ylabel('$\kappa_l$ (kHz)');
    xlim([min_x, max_x]);
    set(gca,'XTick',[]);
    hold off;
    
    subp4 = subplot(4, 2, 7); % form factor
    set(subp4, 'Position', get(subp4, 'Position') .* pos_adj);
    hold on;
    plot(1:n_OPs, intake.C_010, '.', 'Color', SD.black, 'MarkerSize', 20);
    xlabel('operating point');
    ylabel('$C_{010}$ (kHz)');
    xlim([min_x, max_x]);
    hold off;
    
    % receiver paramters:
    subp5 = subplot(4, 2, 2);
    set(subp5, 'Position', get(subp5, 'Position') .* pos_adj);
    % title('Receiver Intake');
    hold on;
    plot(1:n_OPs, pow2db(intake.G_AMP_1q_pow), '.', 'Color', SD.mygreen, 'MarkerSize', 20);
    ylabel('$G_\mathrm{AMP,1quad}$ (dB)');
    xlim([min_x, max_x]);
    set(gca,'XTick',[]);
    hold off;
    
    subp6 = subplot(4, 2, 4);
    set(subp6, 'Position', get(subp6, 'Position') .* pos_adj);
    hold on;
    plot(1:n_OPs, intake.BW_AMP_MHz, '.', 'Color', SD.mygreen, 'MarkerSize', 20);
    ylabel('$B_\mathrm{AMP}$ (MHz)');
    xlim([min_x, max_x]);
    set(gca,'XTick',[]);
    hold off;
    
    subp7 = subplot(4, 2, 6);
    set(subp7, 'Position', get(subp7, 'Position') .* pos_adj);
    hold on;
    plot(1:n_OPs, intake.BW_AMP_MHz .* sqrt(intake.G_AMP_1q_pow), '.', 'Color', SD.mygreen, ...
        'MarkerSize', 20);  % amplitude gain-bandwidth product
    ylabel('$\sqrt{G}BP$ (MHz)');
    xlim([min_x, max_x]);
    set(gca,'XTick',[]);
    hold off;
    
    subp8 = subplot(4, 2, 8);
    set(subp8, 'Position', get(subp8, 'Position') .* pos_adj);
    hold on;
    plot(1:n_OPs, intake.T_fridge_K*1E3, '.', 'Color', SD.mypurple, 'MarkerSize', 20);
    xlabel('operating point');
    ylabel('$T_f$ (mK)');
    xlim([min_x, max_x]);
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Intake Parameters.jpg'], 'jpeg');
        disp('intake parameters figure saved successfully');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate and plot cavity tuning and drift
    
    cav_tune_tx1_kHz = 1E6 * diff(intake.f_cav_tx1_GHz);
    cav_tune_tx2_kHz = 1E6 * diff(intake.f_cav_tx2_GHz);
    OP_half_inds = (1:n_OPs-1) + 0.5;
    fac = 1.1;
    min_y_tune = fac * min([cav_tune_tx1_kHz; cav_tune_tx2_kHz]);
    if min_y_tune > 0
        min_y_tune = min_y_tune / fac^2;
    end
    max_y_tune = fac * max([cav_tune_tx1_kHz; cav_tune_tx2_kHz]);
    if max_y_tune < 0
        max_y_tune = max_y_tune / fac^2;
    end
    min_y_drift = fac * min(cav_drift_tx_kHz);
    if min_y_drift > 0
        min_y_drift = min_y_drift / fac^2;
    end
    max_y_drift = fac * max(cav_drift_tx_kHz);
    if max_y_drift < 0
        max_y_drift = max_y_drift / fac^2;
    end
    
    
    
    fig = figure;
    
    % cavity drift plot
    subp1 = subplot(2, 1, 1);
    % title('Cavity Tuning and Drift');
    %set(subp1, 'Position', get(subp1, 'Position') .* pos_adj);
    hold on;
    p1 = plot(OP_half_inds, cav_tune_tx1_kHz, '.', 'Color', SD.myblue, 'MarkerSize', 35);
    p1a = plot(OP_half_inds, cav_tune_tx1_kHz, '-', 'Color', SD.myblue);
    p2 = plot(OP_half_inds, cav_tune_tx2_kHz, '.', 'Color', SD.myred, 'MarkerSize', 35);
    p2a = plot(OP_half_inds, cav_tune_tx2_kHz, '-', 'Color', SD.myred);
    p3 = plot([min_x, max_x], mean([cav_tune_tx1_kHz; cav_tune_tx2_kHz])*[1,1], '--', 'Color', SD.mypink, ...
        'LineWidth', 1.25);
    p4 = plot([min_x, max_x], median([cav_tune_tx1_kHz; cav_tune_tx2_kHz])*[1,1], '--', 'Color', SD.mypurple, ...
        'LineWidth', 1.25);
    ylabel('$f_\mathrm{step}$ (kHz)');
    xlim([min_x, max_x]);
    ylim([min_y_tune, max_y_tune]);
    xticks([]);
    mylegend([p1, p2, p3, p4], {'pre-', 'post-', 'mean', 'median'}, 'Location', 'Best');
    hold off;
    
    % cavity drift plot
    subp2 = subplot(2, 1, 2);
    %set(subp2, 'Position', get(subp2, 'Position') .* pos_adj);
    hold on;
    p1 = plot(1:n_OPs, cav_drift_tx_kHz, '.', 'Color', SD.myorange, 'MarkerSize', 35);
    p1a = plot(1:n_OPs, cav_drift_tx_kHz, '-', 'Color', SD.myorange, 'LineWidth', 2);
    p2 = plot([min_x, max_x], mean(cav_drift_tx_kHz)*[1,1], '--', 'Color', SD.mypink, ...
        'LineWidth', 1.25);
    p3 = plot([min_x, max_x], median(cav_drift_tx_kHz)*[1,1], '--', 'Color', SD.mypurple, ...
        'LineWidth', 1.25);
    xlabel('operating point');
    ylabel('$f_\mathrm{cav}^\mathrm{post} - f_\mathrm{cav}^\mathrm{pre}$ (kHz)');
    xlim([min_x, max_x]);
    ylim([min_y_drift, max_y_drift]);
    mylegend([p1, p2, p3], {'data', 'mean', 'median'}, 'Location', 'Best');
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Cavity Tuning and Drift.jpg'], 'jpeg');
        disp('cavity tuning and drift figure saved successfully');
    end
end

%% ============================================================================================== %%
%% CHOP TAILS
% the data that falls outside of the analysis band can be sacrificed. However, a buffer will be kept
% on the outer end of the data, where it can be kept, for better SG filtering
IF_cuts.raw_buff_unnorm_IF = raw_IF_mult(:, 1:n_processing_p_buff_bins_IF);
f_IF_buff_MHz = intake.f_IF_MHz(intake.f_IF_MHz <= max_processing_freq_p_buff_IF_MHz);
f_IF_MHz = f_IF_buff_MHz(intake.f_IF_MHz <= max_processing_freq_IF_MHz);
f_IF_min_MHz = min(f_IF_MHz); 
f_IF_max_MHz = max(f_IF_MHz); 

%% ============================================================================================== %%
%% SCALE RAW SPECTRA
% take the RF spectra whose tails have been chopped off and scale them to all have the same mean, as
% in ref. [19]
raw_means = nanmean(IF_cuts.raw_buff_unnorm_IF, 2);
IF_cuts.raw_buff_IF = IF_cuts.raw_buff_unnorm_IF ./ raw_means;

%% ============================================================================================== %%
%% UNFOLD RAW SPECTRA
% generate a frequency vector and an RF spectrum for each "cavity tuning." The frequency vector will
% be relative to the frequency of the cavity, but it will not be known which direction it came from.
% I have concluded in ref. [1] that the best way to do this is to construct frequency vectors of
% twice the length of the IF (which will have a gap in the middle since the IF frequencies do not
% exetend all the way down to DC), and corresponding spectrum vectors. These spectra will be
% symmetric about DC, as that is the best I know how to do with the information the JPA affords us.

% start with an IF vec from 0:
f_IF_buff_from_zero_MHz = 0:raw_res_Hz*1E-6:max(f_IF_buff_MHz);
f_IF_buff_unfold_MHz = [-fliplr(f_IF_buff_from_zero_MHz(2:end)), f_IF_buff_from_zero_MHz];
f_IF_from_zero_MHz = 0:raw_res_Hz*1E-6:f_IF_max_MHz;
f_IF_unfold_MHz = [-fliplr(f_IF_from_zero_MHz(2:end)), f_IF_from_zero_MHz];
% assemble the RF frequency vectors:
f_RF_buff_GHz = f_IF_buff_unfold_MHz*1E-3 + f_LO_GHz;
f_RF_GHz = f_IF_unfold_MHz*1E-3 + f_LO_GHz;
min_f_RF_GHz = min(f_RF_GHz, [], 'all'); % determine min RF frequency
max_f_RF_GHz = max(f_RF_GHz, [], 'all'); % determine max RF frequency
% construct the unfolded spectra with the appropriate gap (full of NaN's):
nan_filler_mat = NaN * ones(n_OPs, n_noisy_center_bins); % will be useful throughout
nan_filler_vec = nan_filler_mat(1, :); % vector version likewise
raw_buff_uncut_RF = [fliplr(IF_cuts.raw_buff_IF), nan_filler_mat, IF_cuts.raw_buff_IF];

% interpolate added noise and squeezing values at all relevant IF frequencies
G_s_IF_pow = interp1(intake.f_abg_Hz*1E-6, intake.G_s_pow.', f_IF_from_zero_MHz, ...
    'linear', 'extrap').';
N_Hp_IF_1q_qta = interp1(intake.f_NH_Hz*1E-6, intake.N_Hp_1q_qta.', f_IF_from_zero_MHz, ...
    'linear', 'extrap').';
G_s_unfold_pow = [fliplr(G_s_IF_pow(:, 2:end)), G_s_IF_pow];
N_Hp_unfold_1q_qta = [fliplr(N_Hp_IF_1q_qta(:, 2:end)), N_Hp_IF_1q_qta];

%% ============================================================================================== %%
%% IF CUTS
% flag and remove IF bins based on whether they exceed a threshold normalized power. Remove
% a specified number of neighboring bins to either side of every contiguous set of flagged bins as
% well. Once bins are removed, the standard deviation of normalized power will decrease, so the
% process repeats iteratively. Should take ~2-3 iterations to converge. Do this whole process
% separately for separate batches of spectra

disp(' ');
disp('data processing begun'); % this marks the start of the data processing

for b = 1:n_batches
    % turn off shortcutting if the relevant file doesn't exist:
    if ~exist([save_to_folder, f, save_to_IF_cuts], 'file') && shortcut_IF_cuts
        warning('no file exists - cannot shortcut IF cuts');
        shortcut_IF_cuts = 0;
    end
    
    if ~shortcut_IF_cuts
        
        this_batch = (intake.batch == b);
        IF_cuts_b{b}.raw_buff_unnorm_IF = IF_cuts.raw_buff_unnorm_IF(this_batch,:); %#ok<*SAGROW>
        IF_cuts_b{b}.raw_buff_IF = IF_cuts.raw_buff_IF(this_batch,:);
        
        disp(['performing IF cuts for batch ', num2str(b), '...']);
        IF_cuts_b{b}.flag_inds = NaN;
        IF_cuts_b{b}.flags = false(1, n_processing_p_buff_bins_IF); % init false
        % track some processing info:
        IF_cuts_b{b}.iters = 0; % init 0
        while ~isempty(IF_cuts_b{b}.flag_inds)
            % change flagged bins from last iteration to NaN's:
            if IF_cuts_b{b}.iters > 0 % won't apply first time through while loop
                IF_cuts_b{b}.raw_buff_IF(:, IF_cuts_b{b}.flag_inds) = NaN;
            end
            
            % normalize:
            % take the IF-mean of the raw spectra:
            IF_cuts_b{b}.mean_raw_buff_IF = nanmean(IF_cuts_b{b}.raw_buff_IF);
            % estimate the spectral baseline with a SG filter capable of omiting NaNs:
            warning('disabling warnings for SG filtering');
            ws = warning('off', 'all');
            [IF_cuts_b{b}.mean_raw_filt_buff_IF, IF_cuts_b{b}.mean_raw_filt_buff_nanless_IF] = ...
                FUNC_sgolayfilt_omitnan_v7(IF_cuts_b{b}.mean_raw_buff_IF, SG_raw_filt_deg, ...
                2 * SG_raw_filt_W_bins + 1, freq_domain_SG);
            ws = warning('on', 'all');
            
            % normalize by the spectral baseline:
            IF_cuts_b{b}.mean_norm_buff_IF = IF_cuts_b{b}.mean_raw_buff_IF ./ IF_cuts_b{b}.mean_raw_filt_buff_IF;
            % subtract off 1 to make mean-0 (not subtracting off the sample mean):
            IF_cuts_b{b}.mean_norm_exc_buff_IF = IF_cuts_b{b}.mean_norm_buff_IF - 1;
            
            % on first iteration, store the raw mean:
            if ~IF_cuts_b{b}.iters
                IF_cuts_b{b}.mean_raw_uncut_buff_IF = IF_cuts_b{b}.mean_raw_buff_IF;
                IF_cuts_b{b}.mean_raw_filt_uncut_buff_IF = IF_cuts_b{b}.mean_raw_filt_buff_IF;
            end
            
            % flag bad bins:
            % flag using the st. dev. of the distribution itself, not its expected value (another slight
            % deviation from Ben)
            IF_cuts_b{b}.flags = IF_cuts_b{b}.mean_norm_exc_buff_IF > IF_discard_thresh_sigma * ...
                nanstd(IF_cuts_b{b}.mean_norm_exc_buff_IF);
            % flag sets of bins adjacent to the original perpetraters:
            [IF_cuts_b{b}.flags, IF_cuts_b{b}.flag_inds] = FUNC_flag_adjacent_v3(IF_cuts_b{b}.flags, n_adj_flag_IF);
            % note that I am choosing to remove the surrounding bins within the iterative loop. This has
            % slightly different results than if I choose to do it afterwards. Since in both cases the
            % bins being removed are of marginal interest, I do not think it matters so much. % note
            % also some of these may be redundant with bins flagged on other iterations. The code
            % handles that fine, but the user should be aware that the total number of cut bins is not
            % necessarily the sum of the number flagged on each iteration
            
            % update tallies/running sums:
            IF_cuts_b{b}.iters = IF_cuts_b{b}.iters + 1;
        end
        % after the last iteration, store the uncut normalized mean as the stored raw spectrum divided
        % by the final SG filter. This makes excesses appear more accurately when plotted after all the
        % iterations than in previous versions (this is new as of v15)
        IF_cuts_b{b}.mean_norm_uncut_buff_IF = ...
            IF_cuts_b{b}.mean_raw_uncut_buff_IF ./ IF_cuts_b{b}.mean_raw_filt_buff_nanless_IF;
        IF_cuts_b{b}.mean_norm_exc_uncut_buff_IF = IF_cuts_b{b}.mean_norm_uncut_buff_IF - 1;
        
        IF_cuts_b{b}.flag_inds = find(isnan(IF_cuts_b{b}.mean_raw_buff_IF));
        % the standard deviation here is really just tracking varaition between the spectra, which for
        % very large sample sizes dominates variation within the spectra. I.e. for us, these spectra are
        % no longer being drawn from the same parent distribution
        IF_cuts_b{b}.stdev_raw_spectrum = nanstd(IF_cuts_b{b}.raw_buff_IF);
        IF_cuts_b{b}.n_flags = sum(isnan(IF_cuts_b{b}.mean_raw_buff_IF(1:n_processing_bins_IF)));
        IF_cuts_b{b}.pct_IF_bins_flagged = 100 * IF_cuts_b{b}.n_flags / n_processing_bins_IF;
        
        % unfold the various spectra now that cuts have been made:
        IF_cuts_b{b}.raw_buff_RF = [fliplr(IF_cuts_b{b}.raw_buff_IF), nan_filler_mat(this_batch,:), IF_cuts_b{b}.raw_buff_IF];
        IF_cuts_b{b}.mean_raw_filt_buff_RF = [fliplr(IF_cuts_b{b}.mean_raw_filt_buff_IF), ...
            nan_filler_vec, IF_cuts_b{b}.mean_raw_filt_buff_IF];
        IF_cuts_b{b}.mean_raw_filt_buff_nanless_RF = [fliplr(IF_cuts_b{b}.mean_raw_filt_buff_nanless_IF), ...
            nan_filler_vec, IF_cuts_b{b}.mean_raw_filt_buff_nanless_IF];
        IF_cuts_b{b}.mean_raw_filt_uncut_buff_RF = [fliplr(IF_cuts_b{b}.mean_raw_filt_uncut_buff_IF), ...
            nan_filler_vec, IF_cuts_b{b}.mean_raw_filt_uncut_buff_IF];
        IF_cuts_b{b}.mean_norm_buff_RF = [fliplr(IF_cuts_b{b}.mean_norm_buff_IF), nan_filler_vec, ...
            IF_cuts_b{b}.mean_norm_buff_IF];
        IF_cuts_b{b}.mean_norm_uncut_buff_RF = [fliplr(IF_cuts_b{b}.mean_norm_uncut_buff_IF), ...
            nan_filler_vec, IF_cuts_b{b}.mean_norm_uncut_buff_IF];
        IF_cuts_b{b}.mean_norm_exc_buff_RF = [fliplr(IF_cuts_b{b}.mean_norm_exc_buff_IF), ...
            nan_filler_vec, IF_cuts_b{b}.mean_norm_exc_buff_IF];
        IF_cuts_b{b}.mean_norm_exc_uncut_buff_RF = [fliplr(IF_cuts_b{b}.mean_norm_exc_uncut_buff_IF), ...
            nan_filler_vec, IF_cuts_b{b}.mean_norm_exc_uncut_buff_IF];
        
        if save_IF_cuts_to_wkspc && b == n_batches % to save time on future runs of the code
            disp('saving IF cut info to workspace...');
            save([save_to_folder, f, save_to_IF_cuts], varname(IF_cuts_b));
            disp([tab, 'IF cut info saved']);
        end
        
    else % to save time, load from workspace set up on a previous run of code
        if b == n_batches
            disp('loading IF cut info from workspace...');
            load([save_to_folder, f, save_to_IF_cuts]);
            disp([tab, 'IF cut info loaded']);
        end
    end
    
end

% display cut info:
for b = 1:n_batches
    disp(['batch ', num2str(b), ' results: '])
    disp([tab, num2str(IF_cuts_b{b}.n_flags), '/', num2str(n_processing_p_buff_bins_IF), ...
        ' IF bins flagged and cut (', num2str(round(IF_cuts_b{b}.pct_IF_bins_flagged, 2)), '% of bins)']);
    disp([tab, 'IF spikes above ', num2str(IF_discard_thresh_sigma), ' sigma cut']);
end


zs = zeros(n_OPs, n_processing_p_buff_bins_IF);
zs2 = zeros(n_OPs, n_processing_p_buff_bins_RF);
IF_cuts.raw_buff_IF = zs;
IF_cuts.mean_raw_buff_IF = zs;
IF_cuts.mean_raw_filt_buff_IF = zs;
IF_cuts.mean_norm_buff_IF = zs;
IF_cuts.mean_norm_exc_buff_IF = zs;
IF_cuts.mean_raw_uncut_buff_IF = zs;
IF_cuts.mean_raw_filt_uncut_buff_IF = zs;
IF_cuts.mean_norm_uncut_buff_IF = zs;
IF_cuts.mean_norm_exc_uncut_buff_IF = zs;
IF_cuts.raw_buff_RF = zs2;
IF_cuts.mean_raw_filt_buff_RF = zs2;

for b = 1:n_batches
    this_batch = (intake.batch == b);
    
    IF_cuts.mean_raw_buff_IF(this_batch,:) = ones(sum(this_batch),1)*IF_cuts_b{b}.mean_raw_buff_IF;
    IF_cuts.mean_raw_filt_buff_IF(this_batch,:) = ones(sum(this_batch),1)*IF_cuts_b{b}.mean_raw_filt_buff_IF;
    IF_cuts.mean_norm_buff_IF(this_batch,:) = ones(sum(this_batch),1)*IF_cuts_b{b}.mean_norm_buff_IF;
    IF_cuts.mean_norm_exc_buff_IF(this_batch,:) = ones(sum(this_batch),1)*IF_cuts_b{b}.mean_norm_exc_buff_IF;
    IF_cuts.mean_raw_uncut_buff_IF(this_batch,:) = ones(sum(this_batch),1)*IF_cuts_b{b}.mean_raw_uncut_buff_IF;
    IF_cuts.mean_raw_filt_uncut_buff_IF(this_batch,:) = ones(sum(this_batch),1)*IF_cuts_b{b}.mean_raw_filt_uncut_buff_IF;
    IF_cuts.mean_norm_uncut_buff_IF(this_batch,:) = ones(sum(this_batch),1)*IF_cuts_b{b}.mean_norm_uncut_buff_IF;
    IF_cuts.mean_norm_exc_uncut_buff_IF(this_batch,:) = ones(sum(this_batch),1)*IF_cuts_b{b}.mean_norm_exc_uncut_buff_IF;
    IF_cuts.mean_raw_filt_buff_RF(this_batch,:) = ones(sum(this_batch),1)*IF_cuts_b{b}.mean_raw_filt_buff_RF;
    
    IF_cuts.raw_buff_IF(this_batch,:) = IF_cuts_b{b}.raw_buff_IF;
    IF_cuts.raw_buff_RF(this_batch,:) = IF_cuts_b{b}.raw_buff_RF;
end

%% ============================================================================================== %%
%% DETERMINE SPECTRUM CUTS
% cut whole spectra (OPs) according to the criteria in ref. [25]
batch_tone_heights = NaN(1, n_batches);
batch_tone_height_stds = NaN(1, n_batches);
batch_avg_pow = NaN(1, n_batches);

for i = 1:n_batches % the cuts group naturally by batch
    unbuff_analysis_range = 1:n_analysis_bins_IF(i);
    
    batch_inds = find(intake.batch == i);
    batch_tone_heights(i) = mean(intake.tone_height(batch_inds));
    batch_tone_height_stds(i) = mean(intake.tone_height_std(batch_inds));
    batch_avg_pow(i) = mean(nanmean(IF_cuts.raw_buff_unnorm_IF(batch_inds, unbuff_analysis_range),2));
    sigma_batch_tone_heights(i) = std(intake.tone_height(batch_inds));
    sigma_batch_tone_height_stds(i) = std(intake.tone_height_std(batch_inds));
    sigma_batch_avg_pow(i) = std(nanmean(IF_cuts.raw_buff_unnorm_IF(batch_inds, unbuff_analysis_range),2));
end
% cuts fall into 6 categories: 
n_cut_types = 6; 
% 1. cuts on avg tone heights within a batch (outside of 3 sigma)
indc_tone_avg = 1;
lbl_tone_avg = '$P_\mathrm{tone}$';
% 2. cuts on std tone heights within a batch (outside of 3 sigma)
indc_tone_std = 2;
lbl_tone_std =  '$\sigma_{P_\mathrm{tone}}$';
% 3. cuts on average spectral power within a batch (outside of 3 sigma)
indc_spec_avg = 3;
lbl_spec_avg =  '$P_\mathrm{spec}$';
% 4-5. cuts on cavity drift (outside of some given value - Kelly and I decided on 60 kHz)
cav_drift_cut_kHz = 60;
indc_f_drift_tx = 4; 
lbl_drift_tx = '$\Delta \nu_\mathrm{tx}$';
indc_f_drift_rfl = 5; 
lbl_drift_rfl = '$\Delta \nu_\mathrm{rfl}$';
% 6. measured squeezing > 1 (i.e. amplification, not squeezing)
indc_nosq = 6;  
lbl_nosq = '$\mathrm{noSQ}$';

n_cut_sig = 3; % use 3 sigma for cuts

cuts = NaN(n_cut_types, n_OPs);

for i = 1:n_OPs
    b = intake.batch(i);
    unbuff_analysis_range = 1:n_analysis_bins_IF(b);

    cuts(indc_tone_avg, i) = abs(intake.tone_height(i) - batch_tone_heights(b)) > ...
        (n_cut_sig * sigma_batch_tone_heights(b));
    cuts(indc_tone_std, i) = abs(intake.tone_height_std(i) - batch_tone_height_stds(b)) > ...
        (n_cut_sig * sigma_batch_tone_height_stds(b));
    cuts(indc_spec_avg, i) = ...
        abs(nanmean(IF_cuts.raw_buff_unnorm_IF(i,unbuff_analysis_range)) - batch_avg_pow(b)) > (n_cut_sig * sigma_batch_avg_pow(b));
end
cuts(indc_f_drift_tx,:) = (abs(cav_drift_tx_kHz) > cav_drift_cut_kHz);
cuts(indc_f_drift_rfl,:) = (abs(cav_drift_rfl_kHz) > cav_drift_cut_kHz); 

cuts(indc_nosq, :) = (intake.direct_squeezing_pow < 1); % squeezing < 1 is bad: it is amplification
any_cut = sum(cuts, 1) > 0;
lbl_any_cut = '$\mathrm{any}$';
n_spec_cuts_tot = sum(any_cut);

%% ============================================================================================== %%
%% PLOT SPECTRUM CUTS
if plot_spec_cuts
    fig = figure;
    hold on;
    % display version which addresses the truncation problem with pcolor:
    disp_cuts = [cuts; any_cut];
    disp_cuts = [disp_cuts; NaN(size(disp_cuts(end,:)))];
    disp_cuts = [disp_cuts, NaN(size(disp_cuts(:,end)))];
    pcolor(disp_cuts);
    % title(['Spectrum Cuts (', num2str(n_spec_cuts_tot), ' tot)']);
    xlabel('$\textrm{Spectrum index}$');
    ylabel('$\textrm{Cut type}$'); 
    yticks(0.5:n_cut_types+1.5);
    yticklabels({'',lbl_tone_avg, lbl_tone_std, lbl_spec_avg, lbl_drift_tx, lbl_drift_rfl, ...
        lbl_nosq, lbl_any_cut});
    
    ylim([1,n_cut_types+2]);
    xlim([1,n_OPs+1]);
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Spectrum Cuts.jpg'], 'jpeg');
        disp('spectrum cuts figure saved successfully');
    end
    
end

%% ============================================================================================== %%
%% NORMALIZE AND SUBTRACT
% normalize the raw spectra by the final version of the SG filtered spectrum:
disp('normalizing individual spectra by average baseline...');
SG_filts.norm_buff_IF = IF_cuts.raw_buff_IF ./ IF_cuts.mean_raw_filt_buff_IF;
norm_exc_buff_IF = SG_filts.norm_buff_IF - 1; % subtracted expected mean of 1
norm_exc_buff_RF = [fliplr(norm_exc_buff_IF), nan_filler_mat, norm_exc_buff_IF];

% statistics of the normalized/excess (mean) spectra:
% find the observed and expected statistics of the mean normalized spectrum data:
experimental_mean_norm_spec_mean = nanmean(IF_cuts.mean_norm_buff_IF,'all');
experimental_mean_norm_spec_std = nanmean(nanstd(IF_cuts.mean_norm_buff_IF,[],2));
expected_mean_norm_spec_std = 1 / sqrt(intake.sub_per_raw * n_OPs);
% find the observed and expected statistics of the mean normalized excess data:
experimental_mean_norm_excess_mean = nanmean(IF_cuts.mean_norm_exc_buff_IF,'all');
experimental_mean_norm_excess_std = nanmean(nanstd(IF_cuts.mean_norm_exc_buff_IF,[],2));
expected_mean_norm_excess_std = expected_mean_norm_spec_std; % same as previous
% find the excpected standard deviation of the normalzed excess data
expected_norm_excess_std = expected_mean_norm_excess_std * sqrt(n_OPs);

%% ============================================================================================== %%
%% PLOT ALL RAW IF SPECTRA
% plots all the raw spectra (with cuts removed) together in the IF
if plot_raw_spectra_IF
    fig = figure;
    hold on;
    col_vec_1 = summer(n_OPs_per_batch(1));
    col_vec_2 = winter(n_OPs_per_batch(2));
    col_vec_3 = autumn(n_OPs_per_batch(3));
    
    for i = 1:n_OPs
        switch intake.batch(i)
            case 1
                col = col_vec_1(i,:);
            case 2
                col = col_vec_2(i-n_OPs_per_batch(1),:);
            case 3
                col = col_vec_3(i-sum(n_OPs_per_batch(1:2)),:);
            otherwise
                col = SD.black;
                warning('more batches than expected');
        end 
        plot(f_IF_buff_MHz, IF_cuts.raw_buff_IF(i,:), 'Color', col);
    end
    % title('All Raw Spectra (cuts removed)');
    xlabel('$\nu_\mathrm{IF}$ (MHz)');
    ylabel('${P}_\mathrm{IF}$ (a.u.)');
    
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'All Raw Spectra.jpg'], 'jpeg');
        disp('all raw spectra figure saved successfully');
    end
end

%% ============================================================================================== %%
%% PLOT RAW RF SPECTRUM AND MEAN FILTER
% plot with and without cuts, and plot alongside the average spectrum. Note that the frequency axis,
% being in the RF, does not exactly correspond to this average if we regard our cavity as moving in
% frequency. But I think it is still informative to plot, as it is the SG filter we are dividing by
if plot_raw_spectrum_RF
    batch_no = intake.batch(indiv_OP_ind);
    
    % use the randoly chosen spectrum:
    spec_to_plot_uncut = raw_buff_uncut_RF(indiv_OP_ind, :); % version without IF cuts
    spec_to_plot = IF_cuts.raw_buff_RF(indiv_OP_ind, :); % overlay version with IF cuts made
    freqs_to_plot = f_RF_buff_GHz(indiv_OP_ind, :); % RF frequencies of this spectrum
    % set axes:
    min_y = 0.9 * min([spec_to_plot_uncut, IF_cuts_b{batch_no}.mean_raw_filt_buff_RF]);
    max_y = 1.01 * max(abs([spec_to_plot_uncut, IF_cuts_b{batch_no}.mean_raw_filt_buff_RF]));
    ax = [min(freqs_to_plot), max(freqs_to_plot), min_y, max_y];
    
    fig = figure;
    hold on;
    %p1 = plot(freqs_to_plot, spec_to_plot_uncut, 'LineWidth', 1.5, 'Color', SD.red);
    p2 = plot(freqs_to_plot, spec_to_plot, 'LineWidth', 1.5, 'Color', SD.black);
    p3 = plot(freqs_to_plot, IF_cuts_b{batch_no}.mean_raw_filt_buff_RF, 'LineWidth', 1.5, ...
        'Color', SD.myblue);
%     lgd_ord = [p2, p1, p3];
%     lgd_labels = {'raw', ['IF cuts (', num2str(IF_cuts_b{batch_no}.n_flags) ,' bins)'], 'filtered mean'};
    lgd_ord = [p2, p3];
    lgd_labels = {'raw spectrum', 'SG filter'};    


    axis(ax);
    % title(['Raw RF Spectrum No. ', num2str(indiv_OP_ind)]);
    xlabel('$\nu_\mathrm{RF}$ (GHz)');
    ylabel('${P}_\mathrm{RF}$ (a.u.)');
    mylegend(lgd_ord, lgd_labels);
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Individual Raw Spectrum.jpg'], 'jpeg');
        disp('individual raw spectrum figure saved successfully');
    end
end

%% ============================================================================================== %%
%% PLOT MEAN IF SPECTRA
if plot_various_mean_spectra_IF
    for b = 1:n_batches 
        
        min_f_IF_buff_MHz = min(f_IF_buff_MHz);
        max_f_IF_anal_MHz = max(max_analysis_freq_IF_MHz(b));
        [~, ind_f_min] = min(abs(min_f_IF_buff_MHz - f_IF_buff_MHz));
        [~, ind_f_max] = min(abs(max_f_IF_anal_MHz - f_IF_buff_MHz));
        range_f = ind_f_min:ind_f_max;
        
        % plot the mean raw spectrum with (foreground) and without (background) the bad IF bins removed:
        f_IF_MHz_omitnan = f_IF_buff_MHz;
        f_IF_MHz_omitnan(IF_cuts_b{b}.flag_inds) = []; % remove the NaN's
        mean_raw_omitnan = IF_cuts_b{b}.mean_raw_buff_IF;
        mean_raw_omitnan(IF_cuts_b{b}.flag_inds) = []; % remove the NaN's
        stdev_raw_omitnan = IF_cuts_b{b}.stdev_raw_spectrum;
        stdev_raw_omitnan(IF_cuts_b{b}.flag_inds) = []; % remove the NaN's
        raw_omitnan_m1sig = mean_raw_omitnan - stdev_raw_omitnan;
        raw_omitnan_p1sig = mean_raw_omitnan + stdev_raw_omitnan;
        ax = [min_f_IF_buff_MHz, max_f_IF_anal_MHz, ...
            0.94 * min(raw_omitnan_m1sig(range_f)), 1.02 * max(raw_omitnan_p1sig(range_f))];
        % region to fill with
        X_reg = [f_IF_MHz_omitnan, fliplr(f_IF_MHz_omitnan)];
        Y_reg = [raw_omitnan_m1sig, fliplr(raw_omitnan_p1sig)];
        % make a version of the no_IF_cut traces that only has data for where the cuts are, and NaN's
        % elsewhere. This is because, the way the smoothing works, the SD.red shows through on the plot
        % if this is not used
        mean_raw_uncut_filt_invertnan = IF_cuts_b{b}.mean_raw_filt_uncut_buff_IF;
        mean_raw_uncut_filt_invertnan(setdiff(1:end, IF_cuts_b{b}.flag_inds)) = NaN;
        
        fig = figure;
        hold on;
        p1 = fill(X_reg, Y_reg, 'y');
        p2 = plot(f_IF_MHz_omitnan, raw_omitnan_m1sig, 'LineWidth', 0.75, 'Color', SD.black);
        p3 = plot(f_IF_MHz_omitnan, raw_omitnan_p1sig, 'LineWidth', 0.75, 'Color', SD.black);
        p4 = plot(f_IF_buff_MHz, IF_cuts_b{b}.mean_raw_uncut_buff_IF, 'LineWidth', 1.5, ...
            'Color', SD.red');
        p5 = plot(f_IF_buff_MHz, IF_cuts_b{b}.mean_raw_buff_IF, 'LineWidth', 1.5, 'Color', SD.myblue);
        p6 = plot(f_IF_buff_MHz, mean_raw_uncut_filt_invertnan, 'LineWidth', 2.5, 'Color', SD.red');
        p7 = plot(f_IF_buff_MHz, IF_cuts_b{b}.mean_raw_filt_buff_IF, 'LineWidth', 2.5, ...
            'Color', SD.black);
        axis(ax);
        % title(['Mean IF Spectrum ', num2str(b), '/', num2str(n_batches)])
        xlabel('$\nu_\mathrm{IF}$ (MHz)');
        ylabel('$\bar{P}_\mathrm{IF}$ (a.u.)');
        mylegend([p5, p4, p7, p6, p1], {'mean raw', ['IF cuts (', num2str(IF_cuts_b{b}.n_flags) , ...
            ' bins)'], 'mean filter', 'mean filter (before cuts)', 'raw $\pm 1 \sigma$ range'});
        hold off;
        
        if save_figs
            saveas(fig, [save_to_folder, f, 'Mean Raw IF Spectrum ', num2str(b), '.jpg'], 'jpeg');
            disp(['mean raw IF spectrum ', num2str(b), '/', num2str(n_batches), ' figure saved successfully']);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % similar idea, but for the mean normalized spectrum:
        ax = [min_f_IF_buff_MHz, max_f_IF_anal_MHz, ...
            0.995 * min(IF_cuts_b{b}.mean_norm_uncut_buff_IF(range_f)), ...
            1.005 * max(IF_cuts_b{b}.mean_norm_uncut_buff_IF(range_f))];
        % make a version of the no_IF_cut traces that only has data for where the cuts are, and NaN's
        % elsewhere. This is because, the way the smoothing works, the SD.red shows through on the plot
        % if this is not used
        mean_norm_uncut_invertnan = IF_cuts_b{b}.mean_norm_uncut_buff_IF;
        mean_norm_uncut_invertnan(setdiff(1:end, IF_cuts_b{b}.flag_inds)) = NaN;
        
        % generate a line at the standard deviation and at the expected standard deviation:
        sigma_experimental_line_lower = experimental_mean_norm_spec_mean - ...
            experimental_mean_norm_spec_std * ones(1, n_processing_p_buff_bins_IF);
        sigma_experimental_line_upper = experimental_mean_norm_spec_mean + ...
            experimental_mean_norm_spec_std * ones(1, n_processing_p_buff_bins_IF);
        sigma_expected_line_lower = 1 + expected_mean_norm_spec_std * ...
            ones(1, n_processing_p_buff_bins_IF);
        sigma_expected_line_upper = 1 - expected_mean_norm_spec_std * ...
            ones(1, n_processing_p_buff_bins_IF);
        
        fig = figure;
        hold on;
        p1 = plot(f_IF_buff_MHz, mean_norm_uncut_invertnan, 'LineWidth', 0.75, 'Color', SD.red');
        p2 = plot(f_IF_buff_MHz, IF_cuts_b{b}.mean_norm_buff_IF, 'LineWidth', 0.75, ...
            'Color', SD.myblue);
        p3 = plot(f_IF_buff_MHz, sigma_experimental_line_lower, 'LineWidth', 1.5, 'Color', SD.mygray);
        p4 = plot(f_IF_buff_MHz, sigma_experimental_line_upper, 'LineWidth', 1.5, 'Color', SD.mygray);
        p5 = plot(f_IF_buff_MHz, sigma_expected_line_lower, '--', 'LineWidth', 1.5, 'Color', SD.black);
        p6 = plot(f_IF_buff_MHz, sigma_expected_line_upper, '--', 'LineWidth', 1.5, 'Color', SD.black);
        axis(ax);
        % title(['Normalized Mean IF Spectrum ', num2str(b), '/', num2str(n_batches)])
        xlabel('$\nu_\mathrm{IF}$ (MHz)');
        ylabel('$\bar{P}_\mathrm{IF} / \bar{P}_\mathrm{IF}^\mathrm{filt}$');
        mylegend([p2, p1, p3, p5], {'spectrum', ['IF cuts (', num2str(IF_cuts_b{b}.n_flags), ' bins)'], ...
            '$\pm 1 \sigma_\mathrm{meas}$', ...
            '$\pm (n_\mathrm{spec} n_\mathrm{avg}\tau_\mathrm{aq} \Delta \nu)^{-1/2}$'}, ...
            'Location', 'SouthEast');
        hold off;
        
        if save_figs
            saveas(fig, [save_to_folder, f, 'Mean Normalized IF Spectrum ', num2str(b), '.jpg'], 'jpeg');
            disp(['mean normalized IF spectrum figure ', num2str(b), '/', num2str(n_batches), ' saved successfully']);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % and finally for the normalized excess spectra (i.e. same with expectation value subtracted):
        % generate a line at the standard deviation and at the expected standard deviation:
        sigma_experimental_line_lower = experimental_mean_norm_excess_mean - ...
            experimental_mean_norm_excess_std * ones(1, n_processing_p_buff_bins_IF);
        sigma_experimental_line_upper = experimental_mean_norm_excess_mean + ...
            experimental_mean_norm_excess_std * ones(1, n_processing_p_buff_bins_IF);
        sigma_expected_line = expected_mean_norm_excess_std * ...
            ones(1, n_processing_p_buff_bins_IF);
        
        % generate the threshold line so it can be drawn in:
        thresh_line_val = IF_discard_thresh_sigma * experimental_mean_norm_excess_std;
        IF_cut_thresh_line = thresh_line_val * ones(1, n_processing_p_buff_bins_IF);
        
        min_y = 1.02 * min(IF_cuts_b{b}.mean_norm_exc_buff_IF(range_f));
        max_y = 1.25 * max([IF_cuts_b{b}.mean_norm_exc_buff_IF(range_f), thresh_line_val]);
        axes_excess = [min_f_IF_buff_MHz, max_f_IF_anal_MHz, min_y, max_y];
        
        % hisogram properties
        u_bound_edge = 1.00000001 * max(abs(IF_cuts_b{b}.mean_norm_exc_buff_IF(range_f)));
        l_bound_edge = -u_bound_edge; % symmetric histogram
        bins_per_sigma = 16; % set a fixed number of bins per standard deviation
        % use the expected sigma in determining the number of bins:
        n_bins = round(bins_per_sigma * u_bound_edge / experimental_mean_norm_excess_std);
        n_edges = n_bins + 1; % 1 more edge than bin
        edges = linspace(l_bound_edge, u_bound_edge, n_edges); % locations of edges
        edge_spacing = edges(2) - edges(1); % since edges are uniformly spaced, this works
        bin_centers = linspace(l_bound_edge + edge_spacing / 2, u_bound_edge - edge_spacing / 2, ...
            n_bins);
        counts = histcounts(IF_cuts_b{b}.mean_norm_exc_buff_IF(range_f), edges); % get counts in each bin
        % calculate the expected Gaussian from the pure statistical expectation (assuming no axion)
        expected_gaussian = (edge_spacing * sum(~isnan(IF_cuts_b{b}.mean_norm_exc_buff_IF(range_f))) / ...
            (expected_mean_norm_excess_std * sqrt(2 * pi))) * gaussmf(bin_centers, ...
            [expected_mean_norm_excess_std, 0]);
        min_x_hist = -max(counts) * 0.012;
        max_x_hist = max([counts, expected_gaussian]) * 1.05;
        axes_hist = [min_x_hist, max_x_hist, min_y, max_y];
        
        % make a version of the no_IF_cut traces that only has data for where the cuts are, and NaN's
        % elsewhere. This is because, the way the smoothing works, the SD.red shows through on the plot
        % if this is not used
        mean_norm_uncut_invertnan = IF_cuts_b{b}.mean_norm_exc_uncut_buff_IF;
        mean_norm_uncut_invertnan(setdiff(1:end, IF_cuts_b{b}.flag_inds)) = NaN;
        
        fig = figure('Position', [SD.default_left, SD.default_bot, SD.default_width * 1.25, ...
            SD.default_height]);
        
        pos_subp1 = [0.1430, 0.1430, 0.5020, 0.7539];
        subp1 = subplot('Position', pos_subp1);    hold on;
        p1 = plot(f_IF_buff_MHz, mean_norm_uncut_invertnan, 'LineWidth', 0.75, 'Color', SD.red');
        p2 = plot(f_IF_buff_MHz, IF_cuts_b{b}.mean_norm_exc_buff_IF, 'LineWidth', 0.75, ...
            'Color', SD.myblue);
        p3 = plot(f_IF_buff_MHz, IF_cut_thresh_line, '--', 'LineWidth', 1.5, 'Color', SD.myorange);
        p4 = plot(f_IF_buff_MHz, sigma_experimental_line_lower, 'LineWidth', 1.5, 'Color', SD.mygray);
        p5 = plot(f_IF_buff_MHz, sigma_experimental_line_upper, 'LineWidth', 1.5, 'Color', SD.mygray);
        p6 = plot(f_IF_buff_MHz, sigma_expected_line, '--', 'LineWidth', 1.5, 'Color', SD.black);
        p7 = plot(f_IF_buff_MHz, -sigma_expected_line, '--', 'LineWidth', 1.5, 'Color', SD.black);
        
        axis(axes_excess);
        % title(['Normalized Mean IF Excess ', num2str(b), '/', num2str(n_batches)]);
        xlabel('$\nu_\mathrm{IF}$ (MHz)');
        ylabel(['$\bar{P}_\mathrm{IF} / \bar{P}_\mathrm{IF}^\mathrm{filt} - 1$']);
        mylegend([p2, p1, p3, p4, p6], {'spectrum', ['IF cuts (', num2str(IF_cuts_b{b}.n_flags) , ...
            ' bins)'], ['cut threshold (', num2str(IF_discard_thresh_sigma), '$\sigma$)'], ...
            '$\pm 1 \sigma_\mathrm{meas}$', ...
            '$\pm (n_\mathrm{spec} n_\mathrm{avg}\tau_\mathrm{aq} \Delta \nu)^{-1/2}$'}, ...
            'Location', 'SouthEast');
        hold off;
        
        % include with this plot a histogram of the mean normalized IF excesses:
        subp2 = subplot('Position', [pos_subp1(1) + pos_subp1(3),  pos_subp1(2), ...
            pos_subp1(3) * 0.55, pos_subp1(4)]);
        
        hold on;
        p8 = plot(counts, bin_centers, '.', 'Color', SD.myblue);
        p9 = plot(expected_gaussian, bin_centers, 'Color', SD.black, 'LineWidth', 1);
        p10 = plot([max(counts) / 50, max_x_hist], [thresh_line_val, thresh_line_val], '--', ...
            'LineWidth', 1.5, 'Color', SD.myorange);
        % title('Excess Distribution');
        xlabel('counts');
        mylegend(p9, 'expected $\mathcal{N}$', 'Location', 'SouthEast');
        yticks([]); % no ticks on y-axis
        axis(axes_hist);
        hold off;
        
        if save_figs
            saveas(fig, [save_to_folder, f, 'Mean Normalized IF Excess Spectrum ', num2str(b), '.jpg'], 'jpeg');
            disp(['mean normalized IF excess spectrum figure ', num2str(b), '/', num2str(n_batches), ' saved successfully']);
        end
    end
end

%% ============================================================================================== %%
%% EXCESS HISTOGRAMMING
if plot_hist_mean_excess
    % unfold the excess power spectra into the RF
    norm_exc_RF = [fliplr(norm_exc_buff_IF), nan_filler_mat, norm_exc_buff_IF];
    
    all_excess_data_norm = reshape(norm_exc_RF, [1, n_OPs * ...
        n_processing_p_buff_bins_RF]) / expected_norm_excess_std;
    u_bound_edge = 1.00000001 * max(abs(all_excess_data_norm)); % just above max abs excess
    l_bound_edge = -u_bound_edge; % symmetric histogram
    bins_per_sigma = 3; % set a fixed number of bins per standard deviation
    % use the expected st. dev. in determining the number of bins:
    n_bins = round(bins_per_sigma * u_bound_edge);
    n_edges = n_bins + 1; % 1 more edge than bin
    edges = linspace(l_bound_edge, u_bound_edge, n_edges); % locations of edges
    edge_spacing = edges(2) - edges(1); % since edges are uniformly spaced, this works
    bin_centers = linspace(l_bound_edge + edge_spacing / 2, u_bound_edge - edge_spacing / 2, ...
        n_bins);
    counts = histcounts(all_excess_data_norm, edges); % get counts in each bin
    % calculate the expected distributions from the pure statistics (assuming no axion):
    expected_gaussian = (edge_spacing * sum(~isnan(all_excess_data_norm)) / sqrt(2 * pi)) * ...
        gaussmf(bin_centers, [1, 0]);
    
    % the plot will be logarithmic in counts
    % set up the axes intelligently:
    max_x = u_bound_edge * 1.05;
    max_y = log10(max([counts, expected_gaussian])) * 1.06;
    ax = [-max_x, max_x, 0, max_y];
    
    fig = figure;
    hold on;
    
    p1 = plot(bin_centers, log10(counts), '.', 'Color', SD.myblue, 'MarkerSize', 22);
    lgd_ord = [p1];
    lgd_labels = {'data'};
    % optionally display expected distribution(s):
    if show_precomb_expected_Gaussians
        p2 = plot(bin_centers, log10(expected_gaussian), 'Color', SD.black, 'LineWidth', 1.75);
        % add appropriate entry to legend:
        lgd_ord = [lgd_ord, p2];
        lgd_labels = [lgd_labels, 'expected $\mathcal{N}$'];
    end
    axis(ax);
    % title('Normalized Excess Power Distribution');
    xlabel('normalized power excess');
    ylabel('counts');
    yticks(0:7);
    yticklabels({'$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$', '$10^5$', '$10^6$', '$10^7$'});
    mylegend(lgd_ord, lgd_labels, 'Location', 'NorthWest');
    text(-0.9, 2, ['$\mu$ = ', num2str(round(nanmean(all_excess_data_norm), 4)), newline, ...
        '$\sigma$ = ', num2str(round(nanstd(all_excess_data_norm), 4))], 'Color', SD.myblue);
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Excess Spectra Normalized Histogram.jpg'], 'jpeg');
        disp('excess spectra normalized histogram figure saved successfully');
    end
end

%% ============================================================================================== %%
%% OBTAIN PROCESSED SPECTRA
% turn off shortcutting if the relevant file doesn't exist:
if ~exist([save_to_folder, f, save_to_SG_filters], 'file') && shortcut_SG_filtering
    warning('no file exists - cannot shortcut SG filtering');
    shortcut_SG_filtering = 0;
end

% at some point I have to stop processing the uncut data in parallel with the cut data. This
% seems like a good point to do that
SG_filts.norm_filters_buff_IF = zeros(n_OPs, n_processing_p_buff_bins_IF);
if ~shortcut_SG_filtering
    % remove the structure from each one by dividing out its own SG filter:
    disp('normalizing individual spectra to their own baselines...');
    warning('disabling warnings for SG filtering');
    ws = warning('off', 'all');
    wb = waitbar(0, 'calculating Savitsky-Golay filters...');
    for i = 1:n_OPs
        [SG_filts.norm_filters_buff_IF(i, :), ~] = ...
            FUNC_sgolayfilt_omitnan_v7(SG_filts.norm_buff_IF(i, :), SG_norm_filt_deg, ...
            2 * SG_norm_filt_W_bins + 1, freq_domain_SG);
        waitbar(i/n_OPs, wb);
    end
    close(wb);
    ws = warning('on', 'all');
    disp([tab, 'all ', num2str(n_OPs), ' Savitsky-Golay filters calculated']);
    
    % divide out the SG filters from the spectra:
    SG_filts.norm_filtd_buff_IF = SG_filts.norm_buff_IF ./ SG_filts.norm_filters_buff_IF;
    disp([tab, 'Savitsky-Golay filters divided out from the spectra']);
    
    % obtain the processed spectra by subtracting off 1 (not the sample mean, notably)
    SG_filts.proc_buff_IF = SG_filts.norm_filtd_buff_IF - 1;
    
    % unfold into RF:
    SG_filts.norm_filters_buff_RF = [fliplr(SG_filts.norm_filters_buff_IF), nan_filler_mat, ...
        SG_filts.norm_filters_buff_IF];
    SG_filts.processed_buff_RF = [fliplr(SG_filts.proc_buff_IF), nan_filler_mat, ...
        SG_filts.proc_buff_IF];
    
    if save_SG_filters_to_wkspc % to save time on future runs of the code
        disp('saving Savitsky-Golay filters and processed spectra to workspace...');
        % create directory if necessary:
        if ~exist(save_to_folder, 'dir')
            mkdir(save_to_folder);
        end
        save([save_to_folder, f, save_to_SG_filters], varname(SG_filts));
        disp([tab, 'Savitsky-Golay filters and processed spectra saved']);
    end
    
else % to save time, load from workspace set up on a previous run of code
    disp('loading Savitsky-Golay filters and processed spectra from workspace...');
    load([save_to_folder, f, save_to_SG_filters]);
    disp([tab, 'Savitsky-Golay filters and processed spectra loaded']);
end

%% ============================================================================================== %%
%% PLOT NORMALIZED EXCESS
% plot the selected operating point's spectrum's normalized excess power along with its SG filter
if plot_normalized_excess_RF
    
    % use the randoly chosen spectrum:    
    spec_to_plot = norm_exc_buff_RF(indiv_OP_ind, :);
    filt_to_plot = SG_filts.norm_filters_buff_RF(indiv_OP_ind, :);
    freqs_to_plot = f_RF_buff_GHz(indiv_OP_ind, :); % RF frequencies of this spectrum
    % set axes:
    max_y = 1.01 * max(abs(spec_to_plot));
    min_y = -max_y;
    ax = [min(freqs_to_plot), max(freqs_to_plot), min_y, max_y];
    
    fig = figure;
    hold on;
    p1 = plot(freqs_to_plot, spec_to_plot, 'LineWidth', 1.25, 'Color', SD.black);
    % p2 = plot([axes(1), axes(2)], nanmean(spec_to_plot) * [1, 1], 'LineWidth', 1, 'Color', SD.mypink);
    % p3 = plot([axes(1), axes(2)], nanmean(spec_to_plot) * [0, 0], 'Linewidth', 1, 'Color', SD.mygray);
    p4 = plot(freqs_to_plot, filt_to_plot - 1, 'Color', SD.red);
    
    % lgd_ord = [p1, p4, p2, p3];
    % lgd_labels = {'excess', 'SG filter - 1', '$\mu_\mathrm{meas}$', '0'};
    lgd_ord = [p1, p4];
    lgd_labels = {'raw excess', 'SG filter - 1'};
    
    
    axis(ax);
    % title(['Normalized Power Excess No. ', num2str(indiv_OP_ind)]);
    xlabel('$\nu_\mathrm{RF}$ (GHz)');
    ylabel('$P_\mathrm{RF} / \bar{P}_\mathrm{IF}^\mathrm{filt} - 1$');
    %mylegend(lgd_ord, lgd_labels);
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Individual Normalized Excess.jpg'], 'jpeg');
        disp('individual normalized excess figure saved successfully');
    end
end

%% ============================================================================================== %%
%% UNBUFFER AND QUANTIFY
% remove the buffer that was kept on the high-IF frequency end of the data to this point (for the SG
% filtering) and then quantify the statistics of the (unbuffered) processed spectra
proc_IF = SG_filts.proc_buff_IF(:, 1:n_processing_bins_IF);
% now remove extra from batches as needed if their analysis bands were shorter than the processing
% band, and also from the beginnings: 
[~, ind_start_analysis_band] = min(abs(min_analysis_freq_IF_MHz - f_IF_MHz));
proc_IF(:, 1:ind_start_analysis_band-1) = NaN; 
for i = 1:n_OPs
    b = intake.batch(i);
    % This assumes that extra bins have already been held onto - i.e. that there are ample
    % bins that have been processed to chop off the number we use for buffering. It is okay
    % if there are some more.
    proc_IF(i, n_analysis_bins_IF(b)+1:end) = NaN;
end

% the frequnecy-averaged mean and standard deviations can be compared to their expected values:
proc_f_avgd_mean = nanmean(nanmean(proc_IF)); % first average across spectra, then freq
proc_f_avgd_var = nanmean(nanvar(proc_IF)); % will give an unbiased estiamte of variance
proc_f_avgd_std = sqrt(proc_f_avgd_var); % sqrt to get st. dev.
proc_expected_f_avgd_std = 1 / sqrt(intake.sub_per_raw);

% summarize:
disp([tab, 'processed spctra statistics: frequency-averaged mean = ', num2str(proc_f_avgd_mean), ...
    ', frequency-averaged st. dev. = ', num2str(proc_f_avgd_std)]);
disp([tab, '(expected: mean = 0, st. dev. = ', num2str(proc_expected_f_avgd_std), ...
    ' - not accounting for effects of filtering)']);

%% ============================================================================================== %%
%% MANUAL RF CUTS
% cut segments of some or all spectra in the rf that are clearly problematic. These have undergone
% manual inspection. Tte idea is to perform these before the SG-filtering, as they could strongly
% affect things there (perhaps more than even truncated data is likely to)

% first make a vector of RF frequencies that are one-sided - i.e. that will allow us to make cuts in
% the IF that will be unfolded onto both sides of the spectra
f_RF_side1_GHz = f_LO_GHz + f_IF_MHz*1E-3;
f_RF_side2_GHz = f_LO_GHz - f_IF_MHz*1E-3;

disp('performing manual rf cuts...');
if N_man > 0
    for i = 1:N_man
        f_cut_GHz = f_man_rf_cut_GHz(i,:); % start and end frequency (inclusive)
        cut_OP_inds = man_cut_indices(i,:); % start and end spectrum indices that cut applies to
        man_cuts1 = zeros(size(f_RF_side1_GHz)); % init zeros
        man_cuts2 = zeros(size(f_RF_side1_GHz)); % init zeros
        for j = cut_OP_inds(1):cut_OP_inds(2)
            man_cuts1(j,:) = f_RF_side1_GHz(j,:) >= f_cut_GHz(1) & ...
                f_RF_side1_GHz(j,:) <= f_cut_GHz(2);
            man_cuts2(j,:) = f_RF_side2_GHz(j,:) >= f_cut_GHz(1) & ...
                f_RF_side2_GHz(j,:) <= f_cut_GHz(2);
        end
        man_cuts1 = boolean(man_cuts1);
        man_cuts2 = boolean(man_cuts2);
        
        % perform manual rf cuts:
        proc_IF(man_cuts1 | man_cuts2) = NaN;
    end
    disp([tab, num2str(N_man), ' manual rf cut(s) performed']);
else
    disp([tab, 'no manual rf cuts specified']);
end

%% ============================================================================================== %%
%% AUTOMATIC RF CUTS
% we have to be more careful with cuts in the rf, as rf spikes will look like axions. We will cut
% from individual spectra, and plot all processed spectra with cuts visualized so that all rf cuts
% can be manually inspected.
disp('performing rf cuts...');
RF_cuts.flags = proc_IF / proc_expected_f_avgd_std > RF_discard_thresh_sigma;
RF_cuts.n_init = sum(sum(RF_cuts.flags));

% now flag all bins close to the original flags:
for i = 1:n_OPs
    [RF_cuts.flags(i,:), ~] = FUNC_flag_adjacent_v3(RF_cuts.flags(i,:), n_adj_flag_RF);
end
RF_cuts.n_final = sum(sum(RF_cuts.flags));
disp([tab, num2str(RF_cuts.n_init), ' rf cuts initially identified; ', num2str(RF_cuts.n_final), ...
    ' after flagging adjacent bins']);

% record how many NaNs there were (from IF cuts) in the processed spectra before applying rf cuts.
% Note this can be a smaller number that from the IF cuts * number of operating points because the
% IF cuts number counts those in the buffer
RF_cuts.n_pre_NaNs = sum(sum(isnan(proc_IF)));

% store copies (IF and RF) of the original, un-rf-cut processed spectra
RF_cuts.proc_unRFcut_IF = proc_IF;
RF_cuts.proc_unRFcut_RF = [fliplr(proc_IF), nan_filler_mat, proc_IF];

% perform RF cuts:
proc_IF(RF_cuts.flags) = NaN;
proc_RF = [fliplr(proc_IF), nan_filler_mat, proc_IF];

% record number of NaNs after applying rf cuts:
RF_cuts.n_post_NaNs_by_OP = sum(isnan(proc_IF), 2); % counts manual cuts
RF_cuts.n_post_NaNs = sum(RF_cuts.n_post_NaNs_by_OP);

disp([tab, 'rf spikes above ', num2str(RF_discard_thresh_sigma), ' sigma cut']);

%% ============================================================================================== %%
%% PLOT CUT NUMBERS
% display the number of RF + IF cuts at each operating point
if plot_n_cuts
    col = SD.red;
    x_min = 0;
    x_max = n_OPs + 1;
    
    y_min = 0;
    y_max = max(1, max(RF_cuts.n_post_NaNs_by_OP)) * 1.05;
    
    fig = figure;
    hold on;
    p1 = plot(1:n_OPs, RF_cuts.n_post_NaNs_by_OP, '.', 'MarkerSize', 15, 'Color', col);
    for i = 1:n_OPs
        plot([i,i], [0,RF_cuts.n_post_NaNs_by_OP(i)], '-', 'LineWidth', 1, 'Color', col);
    end
    % line for IF cuts:
    for b = 1:n_batches
        this_batch = find(intake.batch == b);
        % FLAG - this is not displaying right 
        p2 = plot(this_batch, IF_cuts_b{b}.n_flags*ones(size(this_batch)), '-', 'Color', SD.black);
    end
    lgd_ord = [p1, p2];
    lgd_labels = {'$n_\mathrm{IF}+n_\mathrm{rf}$', '$n_\mathrm{IF}$'};
    mylegend(lgd_ord, lgd_labels, 'Location', 'Best');
    xlim([x_min, x_max]);
    ylim([y_min, y_max]);
    xlabel('operating point');
    ylabel('number of cuts');
    % title('Cuts by Operating Point');
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Number of Cuts.jpg'], 'jpeg');
        disp('number of cuts figure saved successfully');
    end
end

%% ============================================================================================== %%
%% PLOT PROCESSED SPECTRUM
% plot the selected operating point's processed (normalized and uniformly rescaled) spectrum:
if plot_proc_spectrum_RF
    % use the selected operating point's spectrum and its frequency vector:
    spectrum_to_plot = proc_RF(indiv_OP_ind, :);
    freqs_to_plot_GHz = f_RF_GHz(indiv_OP_ind, :);
    % generate lines at the expected and observed 1 sigma marks:
    sigma_experimental_line_lower = proc_f_avgd_mean - proc_f_avgd_std * ...
        ones(1, n_processing_bins_RF);
    sigma_experimental_line_upper = proc_f_avgd_mean + proc_f_avgd_std * ...
        ones(1, n_processing_bins_RF);
    
    sigma_expected_line = proc_expected_f_avgd_std * ones(1, n_processing_bins_RF);
    % set axes for the processed spectrum vs. frequency plot:
    max_y = 1.05 * max(abs(spectrum_to_plot));
    axes_proc_spec = [min(freqs_to_plot_GHz), max(freqs_to_plot_GHz), -max_y, max_y];
    
    u_bound_edge = 1.00000001 * max(abs(spectrum_to_plot)); % just above max abs power
    l_bound_edge = -u_bound_edge; % symmetric histogram
    bins_per_sigma = 30; % set a fixed number of bins per standard deviation
    % use the expected sigma in determining the number of bins:
    n_bins = round(bins_per_sigma * u_bound_edge / proc_expected_f_avgd_std);
    n_edges = n_bins + 1; % 1 more edge than bin
    edges = linspace(l_bound_edge, u_bound_edge, n_edges); % locations of edges
    edge_spacing = edges(2) - edges(1); % since edges are uniformly spaced, this works
    bin_centers = linspace(l_bound_edge + edge_spacing / 2, u_bound_edge - edge_spacing / 2, ...
        n_bins);
    counts = histcounts(spectrum_to_plot, edges); % get counts in each bin
    % calculate the expected Gaussian from the pure statistical expectation (assuming no axion)
    expected_gaussian = (edge_spacing * sum(~isnan(spectrum_to_plot)) / ...
        (proc_expected_f_avgd_std * sqrt(2 * pi))) * gaussmf(bin_centers, ...
        [proc_expected_f_avgd_std, 0]);
    
    % note: this is not binomial distributed, because there is a duplication of bins that occurs.
    % Thus, 1/2 the number of counts is binomial (near-Poisson, in practice), but we multiply the
    % error there by 2, for a net factor of sqrt(2). For now, I treat the bionomial errors as
    % Poisson (up to that factor-of-2), which is generally valid for our histograms. Also, it makes
    % much more sense to use the expected Gaussian as the error source
    % counts_err = sqrt(2) * sqrt(expected_gaussian);
    counts_err = sqrt(2) * sqrt(counts);
    resids = counts - expected_gaussian;
    
    axes_hist = [-max(counts) * 0.012, max([counts, expected_gaussian]) * 1.05, -max_y, max_y];
    
    fig = figure('Position', [SD.default_left, SD.default_bot, SD.default_width * 1.25, ...
        SD.default_height]);
    
    pos_subp1 = [0.1430, 0.1430, 0.5020, 0.7539];
    subp1 = subplot('Position', pos_subp1);
    hold on;
    p1 = plot(freqs_to_plot_GHz, spectrum_to_plot, 'LineWidth', 0.75, 'Color', SD.myred);
    p2 = plot(freqs_to_plot_GHz, sigma_experimental_line_lower, 'LineWidth', 1.5, ...
        'Color', SD.mypink);
    p3 = plot(freqs_to_plot_GHz, sigma_experimental_line_upper, 'LineWidth', 1.5, ...
        'Color', SD.mypink);
    p4 = plot(freqs_to_plot_GHz, sigma_expected_line, '--', 'LineWidth', 1.5, 'Color', SD.black);
    p5 = plot(freqs_to_plot_GHz, -sigma_expected_line, '--', 'LineWidth', 1.5, 'Color', SD.black);
    lgd_ord = [p1, p2, p4];
    % lgd_labels = {'processed spectrum', '$\pm \langle 1 \sigma_\mathrm{meas}\rangle_f$', ...
    %     '$\pm (n_\mathrm{avg}\tau_\mathrm{aq} \Delta \nu)^{-1/2}$'};
    lgd_labels = {'processed spectrum', '$\pm \langle \sigma_\mathrm{meas}\rangle_\nu$', ...
        '$\pm \sqrt{\tau \Delta \nu}$'};
    
    axis(axes_proc_spec);
    % title(['\ \ \ \ Proc. Spec. No. ', num2str(indiv_OP_ind)]);
    xlabel('$\nu_\mathrm{RF}$ (GHz)');
    ylabel('Normalized power excess');
    mylegend(lgd_ord, lgd_labels, 'Location', 'SouthEast');
    hold off;
    
    subp2 = subplot('Position', [pos_subp1(1) + pos_subp1(3),  pos_subp1(2), ...
        pos_subp1(3) * 0.55, pos_subp1(4)]);
    hold on;
    p7 = errorbar(counts, bin_centers, [], [], counts_err, counts_err, '.', 'Color', SD.myred);
    lgd_ord = [p7];
    lgd_labels = {'data'};
    if show_precomb_expected_Gaussians
        p8 = plot(expected_gaussian, bin_centers, 'Color', SD.black, 'LineWidth', 1);
        % add appropriate entry to legend:
        lgd_ord = [lgd_ord, p8];
        lgd_labels = [lgd_labels, 'expected $\mathcal{N}$'];
    end
    % title('Power Distribtuion');
    xlabel('Counts');
    mylegend(lgd_ord, lgd_labels, 'Location', 'SouthEast');
    yticks([]); % no ticks on y-axis
    axis(axes_hist);
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Individual Processed Spectrum.jpg'], 'jpeg');
        disp('individual processed spectrum figure saved successfully');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcualte chi^2 per degree of freedom and plot the residuals from the Gaussian histogram above
    
    % ignore the ones with no counts. But since even the low-counteres are not-quite Gaussian, this
    % is only so ideal a figure of merit
    chi2 = sum((resids(counts ~= 0) ./ counts_err(counts ~= 0)).^2);
    nDOF = length(resids(counts ~= 0)); % there are free parameters - the Gaussian is not a fit
    chi2_per_DOF = chi2 / nDOF;
    % there's an alternative way we might do these calculations - see ref. [20]
    counts_duplicated = 1; % boolean
    delta = FUNC_GaussTest_v3(counts, expected_gaussian, counts_duplicated);
    
    x_min = min(bin_centers);
    x_max = max(bin_centers);
    y_max = 1.025 * max([abs(resids + counts_err), abs(resids - counts_err)]);
    y_min = -y_max;
    ax = [x_min, x_max, y_min, y_max];
    
    fig = figure;
    hold on;
    plot([x_min, x_max], [0,0], '--')
    errorbar(bin_centers, resids, counts_err, '.', 'Color', SD.myred, 'MarkerSize', 16);
    axis(ax);
    % title(['Proc. Resids No. ', num2str(indiv_OP_ind),', $\chi^2/n_\mathrm{DOF} = ', ...
        %num2str(chi2_per_DOF, 3), ', \delta =', num2str(delta, 3), '$'])
    ylabel('counts $-$ expected');
    xlabel('normalized power excess');
    hold off;
    if save_figs
        saveas(fig, [save_to_folder, f, 'Processed Spectrum Residuals.jpg'], 'jpeg');
        disp('processed spectrum residuals figure saved successfully');
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now plot the normalized residuals - better for looking for structure
    fig = figure;
    hold on;
    plot([x_min, x_max], [0,0], '--')
    plot(bin_centers, resids./counts_err, '.', 'Color', SD.myred, ...
        'MarkerSize', 16);
    % title(['Proc. Norm. Resids No. ', num2str(indiv_OP_ind)])
    ylabel('(counts $-$ expected)/$\sqrt{\textrm{expected}}$');
    xlabel('normalized power excess');
    xlim([x_min, x_max]);
    hold off;
    if save_figs
        saveas(fig, [save_to_folder, f, 'Processed Spectrum Normalized Residuals.jpg'], 'jpeg');
        disp('processed spectrum normalized residuals figure saved successfully');
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % quantile-quantile plot
    % might as well do this on the IF spectra
    spec_of_int = proc_IF(indiv_OP_ind, :);
    spec_of_int = spec_of_int(~isnan(spec_of_int));
    n_bins_of_int = length(spec_of_int);
    % the Gaussian with the emperically determined standard deviation would have non-unit slope
    % note we are using the st. dev. for the spectum of interest
    emperical_gauss_slope = std(spec_of_int) / proc_expected_f_avgd_std;
    quants = sort(spec_of_int);
    theo_quants = norminv(((1:n_bins_of_int) - 0.5)/n_bins_of_int) * proc_expected_f_avgd_std;
    max_xy = max(abs([theo_quants, quants])) * 1.05;
    x_emp_line = max_xy*[-1, 1];
    y_emp_line = emperical_gauss_slope*x_emp_line;
    
    fig = figure;
    hold on;
    
    
    plot(theo_quants, quants, '.', 'Color', SD.myred, 'MarkerSize', 12);
    plot(max_xy*[-1, 1], max_xy*[-1, 1], '--', 'LineWidth', 1);
    plot(x_emp_line, y_emp_line, '--', 'LineWidth', 1, 'Color', SD.myblue);
    
    % title(['Q-Q No. ', num2str(indiv_OP_ind)]);
    xlabel('expected $\mathcal{N}$ quantiles');
    ylabel('processed quantiles');
    xlim(max_xy*[-1, 1]);
    ylim(max_xy*[-1, 1]);
    xticks(yticks);
    hold off;
    if save_figs
        saveas(fig, [save_to_folder, f, 'Processed Spectrum Q-Q.jpg'], 'jpeg');
        disp('processed spectrum quantile-quantile figure saved successfully');
    end
    
end

%% ============================================================================================== %%
%% PROCESSED SPECTRA HISTOGRAMMING
if plot_hist_all_proc_spec
    % perform a histogramming of the data for the values in all the processed spectra together:
    % for convenience, put all of the processed spectra along a single row, and normalize to the
    % expected standard deviation:
    all_proc_spec_data_norm = reshape(proc_RF, ...
        [1, n_OPs * n_processing_bins_RF]) / proc_expected_f_avgd_std;
    u_bound_edge = 1.00000001 * max(abs(all_proc_spec_data_norm)); % just above max abs power
    l_bound_edge = -u_bound_edge; % symmetric histogram
    bins_per_sigma = 30; % set a fixed number of bins per standard deviation
    % use the expected st. dev. in determining the number of bins:
    n_bins = round(bins_per_sigma * u_bound_edge);
    n_edges = n_bins + 1; % 1 more edge than bin
    edges = linspace(l_bound_edge, u_bound_edge, n_edges); % locations of edges
    edge_spacing = edges(2) - edges(1); % since edges are uniformly spaced, this works
    bin_centers = linspace(l_bound_edge + edge_spacing / 2, u_bound_edge - edge_spacing / 2, ...
        n_bins);
    counts = histcounts(all_proc_spec_data_norm, edges); % get counts in each bin
    % calculate the expected distributions from the pure statistics (assuming no axion):
    expected_gaussian = (edge_spacing * sum(~isnan(all_proc_spec_data_norm)) / sqrt(2 * pi)) * ...
        gaussmf(bin_centers, [1, 0]);
    
    % note the sqrt(2) is because of duplication - same as for the individual processed spectra
    %counts_err = sqrt(2) * sqrt(expected_gaussian);
    counts_err = sqrt(2) * sqrt(counts);
    resids = counts - expected_gaussian;
    
    % goodness of fit:
    chi2 = sum((resids(counts ~= 0) ./ counts_err(counts ~= 0)).^2);
    nDOF = length(resids(counts ~= 0)); % there are free parameters - the Gaussian is not a fit
    chi2_per_DOF = chi2 / nDOF;
    % alternativve method - see ref. [20]
    counts_duplicated = 1; % boolean
    delta = FUNC_GaussTest_v3(counts, expected_gaussian, counts_duplicated);
    
    % the plot will be logarithmic in counts
    % set up the axes intelligently:
    max_x = u_bound_edge * 1.05;
    max_y = log10(max([counts, expected_gaussian])) * 1.06;
    ax = [-max_x, max_x, -0.5, max_y];
    
    p_err = log10(1 + counts_err ./ counts);
    m_err = log10(1 - counts_err ./ counts);
    
    fig = figure;
    hold on;
    p1 = errorbar(bin_centers, log10(counts), m_err, p_err, '.', 'Color', SD.myred, ...
        'MarkerSize', 22);
    lgd_ord = [p1];
    lgd_labels = {'data'};
    % optionally display expected distribution(s):
    if show_precomb_expected_Gaussians
        p2 = plot(bin_centers, log10(expected_gaussian), 'Color', SD.black, 'LineWidth', 1.75);
        % add appropriate entry to legend:
        lgd_ord = [lgd_ord, p2];
        lgd_labels = [lgd_labels, 'expected $\mathcal{N}$'];
    end
    axis(ax);
    % title(['Proc. Distribution, ', '$\chi^2/n_\mathrm{DOF} = ', ...
    %     num2str(chi2_per_DOF, 3), ', \delta =', num2str(delta, 3), '$']);
    xlabel('Renormalized power excess');
    ylabel('Counts');
    if ceil(max_x) <= 10 % in case max_x is very high, do not do this
        xticks(-ceil(max_x):ceil(max_x));
    end
    yticks(0:7);
    yticklabels({'$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$', '$10^5$', '$10^6$', '$10^7$'});
    mylegend(lgd_ord, lgd_labels, 'Location', 'NorthEast');
    text(-1.25, 2, ['$\mu$ = ', num2str(round(nanmean(all_proc_spec_data_norm), 6)), newline, ...
        '$\sigma$ = ', num2str(round(nanstd(all_proc_spec_data_norm), 8))], 'Color', SD.myred);
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Processed Spectra Normalized Histogram.jpg'], 'jpeg');
        disp('processed spectra normalized histogram figure saved successfully');
    end
end

%% ============================================================================================== %%
%% PLOT PROCESSED SPECTRA WAERFALL
if plot_waterfall_spectra_RF
    fig = figure('Position', ...
        [SD.default_left, SD.default_bot, SD.default_width * 1.25, SD.default_height * 1.25]);
    hold on;
    
    % for the sake of the waterfall plot, renormalize s.t. the st. dev. = 1:
    renorm_proc_RF = proc_RF / proc_expected_f_avgd_std;
    
    x_p_min = min(renorm_proc_RF, [], 'all');
    x_p_max = max(renorm_proc_RF, [], 'all') + RF_cutoff * (n_OPs-1);
    
    for s = 1:n_OPs % loop over all the processed spectra
        wf_add = (s-1) * RF_cutoff;
        b = intake.batch(s);
        
        proc2plot = renorm_proc_RF(s, :) + wf_add; % waterfall them
        freqs2plot_GHz = f_RF_GHz(s, :);
        
        % draw some horizontal lines - improves zoom-ability not to have these lines
%         plot([min_f_RF_GHz, max_f_RF_GHz], wf_add * [1, 1], 'Color', SD.mygray, 'LineWidth', 0.75);
%         plot([min_f_RF_GHz, max_f_RF_GHz], wf_add * [1, 1] + 1, 'Color', SD.mygray, 'LineWidth', 0.5);
%         plot([min_f_RF_GHz, max_f_RF_GHz], wf_add * [1, 1] - 1, 'Color', SD.mygray, 'LineWidth', 0.5);
        
        lw = 1;
        switch b % color different batches
            case 1
                col = SD.myred;
            case 2
                col = SD.mygreen;
            case 3
                col = SD.mypurple;
        end
        if any_cut(s)
            col = SD.black; % color cuts black
        end
        if s == indiv_OP_ind
            col = SD.myred;
        end
        
        plot(freqs2plot_GHz, proc2plot, 'LineWidth', lw, 'Color', col);
        
%         fsz = 8;
%         text(max_f_RF_GHz + (max_f_RF_GHz - min_f_RF_GHz) * 0.02, wf_add, ['$\mu = ', ...
%             num2str(round(nanmean(proc2plot - wf_add), 5), 2), '$'], 'FontSize', fsz);
    end
    
    % title('$\textrm{Processed Spectra}$');
    xlabel('$\nu_\mathrm{RF}\ (\textrm{GHz})$');
    %ylabel('$\textrm{renormalized power excess}$');
    ylabel('$\textrm{Normalized power excess}$');
    
    xlim([min_f_RF_GHz, max_f_RF_GHz + (max_f_RF_GHz - min_f_RF_GHz) * 0.0]);
    ylim([x_p_min, x_p_max]);
    
    tick_mult = 100;
    ytick_vals = 0:tick_mult*RF_cutoff:1E4;
    ytick_labs = strings(size(ytick_vals));
    for i = 1:length(ytick_vals)
        ytick_labs(i) = ['+', num2str(ytick_vals(i))];
    end
    %yticks(ytick_vals);
    %yticklabels(ytick_labs);
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Processed Spectra Waterfall'], 'jpeg');
        disp('Processed spectra waterfall plot saved');
    end
    
    hold off;
end

%% ============================================================================================== %%
%% PROCESS CAVITY SPECTROSCOPY
% because the data being intaken from the measurement has already fitted the cavity before and
% after, this section is only used as a check. It is fine to go with the values from the original
% measurement (see ref. [15])
disp('processing cavity spectroscopy data...');
% define the relevant centered, unit-height Lorentzian function:
cav_lor = @(f, w, c) w.^2 ./ (w.^2 + 4*(f-c).^2); % f is frequency, w is FWHM linewidth, c is center

cav_fit_FWHM_tx1_kHz = NaN(n_OPs, 1);
cav_fit_FWHM_tx2_kHz = NaN(n_OPs, 1);
cav_fit_cent_tx1_kHz = NaN(n_OPs, 1);
cav_fit_cent_tx2_kHz = NaN(n_OPs, 1);
if ~bypass_cav_spec
    wb = waitbar(0, 'processing cavity spectroscopy...');
else
    warning('bypassing cavity spectroscopy');
end

for i = 1:n_OPs % perform fit at each operating point
    if ~bypass_cav_spec || i == indiv_OP_ind
        % perform the fit for both the initial and final ENA sweeps for each operating point:
        Ssw2_tx1_pow = abs(intake.ENA_Ssw_tx1).^2; % ENA trace from before raw acquisition
        Ssw2_tx2_pow = abs(intake.ENA_Ssw_tx2).^2; % from after
        Ssw2_norm_tx1_pow = Ssw2_tx1_pow ./ max(Ssw2_tx1_pow, [], 2); % normalize
        Ssw2_norm_tx2_pow = Ssw2_tx2_pow ./ max(Ssw2_tx2_pow, [], 2);
        
        % common to both fits:
        peak_guess = 1; % the normalization ensures only that the max is 1, peak will be close
        
        % fit to first ENA trace:
        cent_guess_kHz = intake.f_cav_tx1_GHz(i) * 1E6;
        FWHM_guess_kHz = intake.cav_kappa_ext_tx1_kHz(i) + intake.cav_kappa_loss_tx1_kHz(i);
        [cav_fit_result_tx1, ~, ~] = FUNC_noncentered_Lorentzian_fit_v2(...
            intake.ENA_f_tx1_GHz(i,:)*1E6, Ssw2_norm_tx1_pow(i,:), FWHM_guess_kHz, peak_guess, ...
            cent_guess_kHz);
        cav_fit_FWHM_tx1_kHz(i) = cav_fit_result_tx1.w; %the FWHM is the parameter we care about
        cav_fit_cent_tx1_kHz(i) = cav_fit_result_tx1.c;
        
        % fit to second ENA trace:
        cent_guess_kHz = intake.f_cav_tx2_GHz(i) * 1E6;
        FWHM_guess_kHz = intake.cav_kappa_ext_tx2_kHz(i) + intake.cav_kappa_loss_tx2_kHz(i);
        [cav_fit_result_tx2, ~, ~] = FUNC_noncentered_Lorentzian_fit_v2(...
            intake.ENA_f_tx2_GHz(i,:)*1E6, Ssw2_norm_tx2_pow(i,:), FWHM_guess_kHz, peak_guess, ...
            cent_guess_kHz);
        cav_fit_FWHM_tx2_kHz(i) = cav_fit_result_tx2.w; % the FWHM is the parameter we care about
        cav_fit_cent_tx2_kHz(i) = cav_fit_result_tx2.c;
        
        
        if i == indiv_OP_ind && plot_cav_spec % for the selected operating point, save info
            cav_plt.f_ENA_IF_tx1_MHz = intake.ENA_f_tx1_GHz(i,:)*1E3 - cav_fit_cent_tx1_kHz(i)*1E-3;
            % IF frequencies are defined relative to the first transmission measurement:
            cav_plt.f_ENA_IF_tx2_MHz = intake.ENA_f_tx2_GHz(i,:)*1E3 - cav_fit_cent_tx1_kHz(i)*1E-3;
            cav_plt.raw_trans_tx1_dB = pow2db(Ssw2_norm_tx1_pow(i,:));
            cav_plt.raw_trans_tx2_dB = pow2db(Ssw2_norm_tx2_pow(i,:));
            cav_plt.lor_fit_tx1_dB = pow2db(cav_lor(f_IF_unfold_MHz*1E3 + cav_fit_cent_tx1_kHz(i), ...
                cav_fit_FWHM_tx1_kHz(i), cav_fit_cent_tx1_kHz(i)));
            % IF frequencies are defined relative to the first transmission measurement:
            cav_plt.lor_fit_tx2_dB = pow2db(cav_lor(f_IF_unfold_MHz*1E3 + cav_fit_cent_tx1_kHz(i), ...
                cav_fit_FWHM_tx2_kHz(i), cav_fit_cent_tx2_kHz(i)));
        end
        if ~bypass_cav_spec
            waitbar(i/n_OPs, wb);
        end
    end
end
if ~bypass_cav_spec
    close(wb);
end
disp([tab, 'cavity spectroscopy data processed...']);

%% ============================================================================================== %%
%% PLOT CAVITY SPECTROSCOPY
if plot_cav_spec
    % plot the cavity spectroscopy against its fit for the selected operating point - the values
    % from the fit in this analysis are not actually what is used in the rescaling, note. They serve
    % only as checks
    max_x = max([cav_plt.f_ENA_IF_tx1_MHz, cav_plt.f_ENA_IF_tx2_MHz]);
    min_y = min([cav_plt.raw_trans_tx1_dB, cav_plt.lor_fit_tx1_dB, cav_plt.raw_trans_tx1_dB, ...
        cav_plt.lor_fit_tx2_dB]);
    ax = [-max_x, max_x, min_y, 0.5];
    
    fig = figure;
    hold on;
    p1 = plot(cav_plt.f_ENA_IF_tx1_MHz, cav_plt.raw_trans_tx1_dB, 'Color', SD.myred, 'LineWidth', 5);
    p1a = plot(cav_plt.f_ENA_IF_tx2_MHz, cav_plt.raw_trans_tx2_dB, 'Color', SD.black, 'LineWidth', 3);
    
    p2 = plot(f_IF_unfold_MHz, cav_plt.lor_fit_tx1_dB, 'Color', SD.myblue, 'LineWidth', 4);
    % IF frequencies are defined relative to the first transmission measurement
    p2a = plot(f_IF_unfold_MHz, cav_plt.lor_fit_tx2_dB, 'Color', SD.gorse, 'LineWidth', 2);
    
    axis(ax);
    xlabel('$\nu_\mathrm{IF}$ (MHz)');
    ylabel('$S_\mathrm{sw}$ (dB)');
    % title('Cavity Spectroscopy');
    mylegend('tx1 data', 'tx2 data', 'tx1 fit', 'tx2 fit');
    text(0, .75*ax(3), ['$\mathrm{FWHM}_1 = ', ...
        num2str(round(cav_fit_FWHM_tx1_kHz(indiv_OP_ind))), ...
        '$ kHz'], 'Color', SD.myred, 'HorizontalAlignment', 'center');
    text(0, .85*ax(3), ['$\mathrm{FWHM}_2 = ', ...
        num2str(round(cav_fit_FWHM_tx2_kHz(indiv_OP_ind))), ...
        '$ kHz'], 'Color', SD.black, 'HorizontalAlignment', 'center');
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Cavity Spectroscopy.jpg'], 'jpeg');
        disp([tab, 'cavity spectroscopy figure saved successfully']);
    end
end

%% ============================================================================================== %%
%% OBTAIN NOISE PROFILE
% the cavity spectroscopy gives the signal power, up to a constant. The noise must be determined
% through a series of measurements, as described in refs. [5-7]. This section is based on the code
% in ref. [8]
disp('rescaling spectra...');
disp([tab, 'obtaining noise profile...']);

% cav measurement port power reflection as a function of detuning from the cavity:
chi_mm2_vec_unfold_pow = ((cav_kappa_diff_kHz*1E-3).^2 + 4*f_IF_unfold_MHz.^2) ./ ...
    ((cav_kappa_tot_avg_kHz*1E-3).^2 + 4*f_IF_unfold_MHz.^2);
% and at the window where we perform our hot rod/squeezing calibration measurements (ref. [18])
chi_mm2_Wdet_pow = ((cav_kappa_diff_kHz*1E-3).^2 + 4*intake.W_center_MHz^2) ./ ...
    ((cav_kappa_tot_avg_kHz*1E-3).^2 + 4*intake.W_center_MHz^2);

S_f_1q_qta = FUNC_singlequad_S_of_T_f_v1(intake.T_fridge_K, f_cav_avg_GHz); % fridge 1-quad spec den
S_r_1q_qta = FUNC_singlequad_S_of_T_f_v1(intake.T_rod_K, f_cav_avg_GHz); % rod 1-quad spec den
% all powers here will be single-quadrature

% HEMT+ and fridge single-quadrature noise powers over the relevant badnwdith:
P_Hp_W = h_Planck_Js * f_cav_avg_GHz*1E9 .* N_Hp_unfold_1q_qta * intake.W_width_MHz*1E6; % added noise
P_f_W = h_Planck_Js * f_cav_avg_GHz*1E9 .* S_f_1q_qta * intake.W_width_MHz*1E6; % fridge input noise


% calculate noise power spectrum entering and exiting cavity:
S_enter_1q_qta = (1-intake.SQ_cav_eta)*S_f_1q_qta + intake.SQ_cav_eta*G_s_unfold_pow.*S_f_1q_qta;
S_exit_unfold_1q_qta = S_enter_1q_qta .* chi_mm2_vec_unfold_pow + ...
    (1 - chi_mm2_vec_unfold_pow) .* intake.S_c_1q_qta;

% now we let the frequency dependence enter (the cavity reflection's was already included, above)
% AMP gain as a function of detuning from the cavity, bottoms out at 1
% FLAG - is intake.G_AMP_1q_pow the right value for on resonance. Is it an off res value?
G_AMP_unfold_1q_pow = 1 + (intake.G_AMP_1q_pow - 1) ./ ...
    ((2 * f_IF_unfold_MHz ./ intake.BW_AMP_MHz).^2 + 1);
S_preAMP_unfold_1q_qta = S_exit_unfold_1q_qta * intake.cav_AMP_eta + ...
    S_f_1q_qta * (1 - intake.cav_AMP_eta) + N_Hp_unfold_1q_qta ./ G_AMP_unfold_1q_pow;

% a version with no AMP-frequency dependance for comparison:
S_preAMP_no_AMP_dep_unfold_1q_qta = S_exit_unfold_1q_qta * intake.cav_AMP_eta + ...
    S_f_1q_qta .* (1 - intake.cav_AMP_eta) + N_Hp_unfold_1q_qta ./ intake.G_AMP_1q_pow;
% flat version that neglects the hot rod and treats the cavity modes as squeezed (unrealistic):
S_preAMP_flat_unfold_1q_qta = S_enter_1q_qta * intake.cav_AMP_eta + S_f_1q_qta * ...
    (1 - intake.cav_AMP_eta) + N_Hp_unfold_1q_qta ./ intake.G_AMP_1q_pow;

disp([tab, 'noise profile obtained...']);

%% ============================================================================================== %%
%% PLOT NOISE PROFILE
if plot_noise_profile
    % plot the noise power profile for the selected individual operating point
    
    % set axes limits:
    all_trace_vec = [S_preAMP_unfold_1q_qta(indiv_OP_ind,:), ...
        S_preAMP_no_AMP_dep_unfold_1q_qta(indiv_OP_ind,:), ...
        S_preAMP_flat_unfold_1q_qta(indiv_OP_ind,:)];
    max_y = 1.025 * max(all_trace_vec);
    min_y = 0.815 * min(all_trace_vec);
    ax = [-max_processing_freq_IF_MHz, max_processing_freq_IF_MHz, min_y, max_y];
    
    fig = figure;
    hold on;
    p1 = plot(f_IF_unfold_MHz, S_preAMP_unfold_1q_qta(indiv_OP_ind,:), '-', 'Color', SD.mygreen);
    p2 = plot(f_IF_unfold_MHz, S_preAMP_no_AMP_dep_unfold_1q_qta(indiv_OP_ind,:), '--', ...
        'Color', SD.mygreen);
    p3 = plot(f_IF_unfold_MHz, S_preAMP_flat_unfold_1q_qta(indiv_OP_ind,:), '-.', ...
        'Color', SD.mygreen);
    axis(ax);
    % title(['Anchored Haloscope Noise No. ', num2str(indiv_OP_ind)]);
    xlabel('$\nu_\mathrm{IF}$ (MHz)');
    ylabel('$S_\mathrm{pre-AMP}$ (qta)');
    legend([p3, p2, p1], ...
        {'$|\chi_{mm}(f)|^2 = 1$, $G_\mathrm{AMP}(f) = G_\mathrm{AMP}^\mathrm{max}$', ...
        '$G_\mathrm{AMP}(f) = G_\mathrm{AMP}^\mathrm{max}$', 'total system noise'}, ...
        'Location', 'South', 'FontSize', 14, 'LineWidth', 2.5);
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Individual Noise Profile.jpg'], 'jpeg');
        disp([tab, 'individual noise profile figure saved successfully']);
    end
end

%% ============================================================================================== %%
%% OBTAIN RESCALING SNR
% we now have the signal and the noise versus frequency. There are 3 IF frequency-dependant
% contributions: the cavity Lorentzian filtering of the signal, the cavity-Lorentzian filtering of
% the unsqueezed (interal loss) noise, and the AMP bandwidth in the presence of non-zero added noise
disp([tab, 'calculating rescaling SNR for a KSVZ axion...']);

% as shown in ref. [3], the Lorentzian we want maximizes at 1. We will use the smooth Lorentzian
% curve that comes from knowing the bandwidth (from the measurement, not from the processing fit)
% FLAG - assuming Lotentzian centers at 0. Might it be better expressed centering at the mean of the
% 2 transmission measurements
cav_trans_norm_fit_unfold_pow = cav_lor(f_IF_unfold_MHz*1E3, cav_kappa_tot_avg_kHz, 0);

% calculate the signal from the cavity Lorentzian, yielding P_ij in ref. [17]:
% note: it is correct to use kappa in Hz here, not rad/s as U0 contains the factor of 2pi (ref. [2])
% note: U0_J already includes the cav-AMP transmission
% FLAG - is P_sig single-quadrature?
P_sig_tot_unfold_W = U0_J * (1E6*f_LO_GHz ./ cav_kappa_tot_avg_kHz).^2 .* ...
    cav_kappa_ext_avg_kHz*1E3 .* intake.C_010 .* cav_trans_norm_fit_unfold_pow; % P_ij, in ref. [17]
% factor of 2 in denominator is for single-quadrature
S_sig_tot_unfold_1q_qta = P_sig_tot_unfold_W ./ (2 * h_Planck_Js * f_LO_GHz*1E9 * raw_res_Hz);

% obtain the signal-to-noise power ratio:
resc_SNR_unfold = S_sig_tot_unfold_1q_qta ./ S_preAMP_unfold_1q_qta; % unfolded version
resc_SNR = resc_SNR_unfold(:, folded_bin_range);

disp([tab, 'rescaling SNR calculated']);

%% ============================================================================================== %%
%% PLOT RESCALING SNR
if plot_resc_SNR
    % plot the SNR used for rescaling for the selected individual operating point
    
    % unfolded SNRs that do not have the entire frequency dependency included, for comparison:
    % 1. includes the cavity Lorentian profile but not the finite AMP BW:
    rescaling_SNR_no_AMP_det_unfold = S_sig_tot_unfold_1q_qta(indiv_OP_ind,:) ./ ...
        S_preAMP_no_AMP_dep_unfold_1q_qta(indiv_OP_ind,:);
    % 2. includes only the cavity Lorentzian profile (noise is treated as fully squeezed - even the
    % internal loss modes; AMP has infinite BW):
    rescaling_SNR_signal_Lorentzian_only_unfold = S_sig_tot_unfold_1q_qta(indiv_OP_ind,:) ./ ...
        S_preAMP_flat_unfold_1q_qta(indiv_OP_ind,:);
    
    % set axes limits:
    all_trace_vec = [resc_SNR_unfold(indiv_OP_ind,:), rescaling_SNR_no_AMP_det_unfold, ...
        rescaling_SNR_signal_Lorentzian_only_unfold];
    max_y = 1.025 * max(all_trace_vec);
    ax = [-max_processing_freq_IF_MHz, max_processing_freq_IF_MHz, 0, max_y];
    
    fig = figure;
    hold on;
    p2 = plot(f_IF_unfold_MHz, rescaling_SNR_no_AMP_det_unfold, '--', 'Color', SD.black,'LineWidth', 1);
    p3 = plot(f_IF_unfold_MHz, rescaling_SNR_signal_Lorentzian_only_unfold, ':', ...
        'Color', SD.black, 'LineWidth', 2);
    p1 = plot(f_IF_unfold_MHz, resc_SNR_unfold(indiv_OP_ind,:), '-', 'Color', SD.myred, 'LineWidth', 2.5);

    axis(ax);
    % title(['Rescaling SNR No. ', num2str(indiv_OP_ind)]);
    xlabel('$\nu_\mathrm{IF}$ (MHz)');
    ylabel('$\alpha^\mathrm{KSVZ}_{1\textrm{-}\mathrm{bin}}$');
    legend([p3, p2, p1], ...
        {'$|\chi_{ml}(\nu)|^2 = 0$, $G_\mathrm{A}(\nu) = G_\mathrm{A}(\nu_p/2)$', ...
        '$G_\mathrm{A}(\nu) = G_\mathrm{A}(\nu_p/2)$', 'total system SNR'}, ...
        'Location', 'South', 'FontSize', 14, 'LineWidth', 2.5);
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Individual Rescaling SNR.jpg'], 'jpeg');
        disp([tab, 'individual rescaling SNR figure saved successfully']);
    end
end

%% ============================================================================================== %%
%% OBTAIN RESCALED SPECTRA
% rescale by taking into account the IF-dependant sensitivity, with the three frequency-dependant
% contributions of the previous several sections

% obtain the rescaled spectrum in the IF first
resc_IF = proc_IF ./ resc_SNR;
% unfold into RF:
resc_RF = [fliplr(resc_IF), nan_filler_mat, resc_IF];

% calculate expected st. dev.
resc_expected_std_IF = proc_expected_f_avgd_std ./ resc_SNR;
resc_expected_std_IF(isnan(proc_IF)) = NaN; % apply all IF and rf cuts to the st devs as well
resc_expected_std_RF = [fliplr(resc_expected_std_IF), nan_filler_mat, resc_expected_std_IF];

disp('spectra rescaled');

%% ============================================================================================== %%
%% PLOT RESCALED SPECTRUM
if plot_resc_spectrum_RF
    % use the selected operating point's spectrum and its frequency vector:
    freqs_to_plot = f_RF_GHz(indiv_OP_ind, :);
    spectrum_to_plot = resc_RF(indiv_OP_ind, :);
    
    resc_abs_max = max(abs(spectrum_to_plot));
    ax = [min(freqs_to_plot), max(freqs_to_plot), -resc_abs_max, resc_abs_max];
    
    fig = figure;
    hold on;
    plot_batch = intake.batch(indiv_OP_ind);
    switch plot_batch
        case 1
            col = SD.myred;
        case 2
            col = SD.mygreen;
        case 3
            col = SD.myred;% myblue
        otherwise
            error('invalid batch specified');
    end
    
    p1 = plot(freqs_to_plot, spectrum_to_plot, 'Color', col, 'LineWidth', 0.75);
    lgd_ord = [p1];
    lgd_labels = {'spectrum'};
    axis(ax);
    xlabel('$\nu_\mathrm{RF}$ (GHz)');
    ylabel('Rescaled power excess');
    % title(['Rescaled Spectrum No. ', num2str(indiv_OP_ind)]);
    % mylegend(lgd_ord, lgd_labels, 'Location', 'South');
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Individual Rescaled Spectrum.jpg'], 'jpeg');
        disp('individual rescaled spectrum figure saved successfully');
    end
end

%% ============================================================================================== %%
%% PERFORM SPECTUM CUTS
% See sec. on determining spectrum cuts, above
for i = 1:n_OPs
    if any_cut(i) % if any of the spectrum cut criteria were met, make the cut
        resc_IF(i,:) = NaN(size(resc_IF(i,:)));
        resc_RF(i,:) = NaN(size(resc_RF(i,:)));
        resc_expected_std_IF(i,:) = NaN(size(resc_expected_std_IF(i,:)));
        resc_expected_std_RF(i,:) = NaN(size(resc_expected_std_RF(i,:)));
    end
end

%% ============================================================================================== %%
%% OBTAIN COMBINED SPECTRA
% beginning at this point, we work fully in the RF
% the combined spectra takes the ML-weighted sum of all contributing spectra into each of its bins
disp('combining spectra...');

% determine the number of bins in the combined spectrum (allow for non-uniform stepping) and the
% combined spectrum frequency vector:
range_f_RF_GHz = max_f_RF_GHz - min_f_RF_GHz;
n_bins_comb = round(range_f_RF_GHz*1E9/raw_res_Hz + 1);
f_comb_GHz = linspace(min_f_RF_GHz, max_f_RF_GHz, n_bins_comb);

% determine the start and stop RF indecies for each row of the matrices below:
start_freqs_GHz = min(f_RF_GHz, [], 2);
stop_freqs_GHz = max(f_RF_GHz, [], 2);
% rounding ensures we get integers for the indecies:
start_inds = round((start_freqs_GHz - min_f_RF_GHz)*1E9 / raw_res_Hz + 1);
stop_inds = round((stop_freqs_GHz - min_f_RF_GHz)*1E9 / raw_res_Hz + 1);

% form matrices of spectra, variances, inverse variances, and SNRs. They are indexed by ref. [2]'s i
% and k, which mathematically is a summation over j, as in ref. [9]:
resc_comb_RF = NaN(n_OPs, n_bins_comb);
var_resc_comb_RF = NaN(n_OPs, n_bins_comb);
inv_var_resc_comb_RF = NaN(n_OPs, n_bins_comb);
SNR_resc_comb_RF = NaN(n_OPs, n_bins_comb);
% fill in these matrices:
for i = 1:n_OPs
    range = start_inds(i):stop_inds(i); % range of indecies to fill in
    
    resc_comb_RF(i, range) = resc_RF(i, :);
    var_resc_comb_RF(i, range) = resc_expected_std_RF(i, :).^2;
    inv_var_resc_comb_RF(i, range) = 1 ./ var_resc_comb_RF(i, range);
    SNR_resc_comb_RF(i, range) = resc_SNR_unfold(i, :);
end
% replace the locations in the rescaling SNR where there are NaN's in the rescaled-combined spectra
% with NaN's. This will prove important in accurately calculating the expected standard deviations:
SNR_resc_comb_RF(isnan(resc_comb_RF)) = NaN;
% form a matrix of weights with the same dimensions as above as in refs. [2] and [9]:
weight_resc_comb_RF = inv_var_resc_comb_RF ./ nansum(inv_var_resc_comb_RF);
% determine number of contributing raw spectra to combined spectrum bins (cuts do not contribute):
comb_cont_raw_spectra = sum(~isnan(resc_comb_RF));

% form the combined spectrum, which averages all contributing spectra, properly weighted (note that
% later the combined spectrum and its standard deviation will be redefined. I denote the
% difference between the original and redefined versions in the variable names):
comb_pre_sum = weight_resc_comb_RF .* resc_comb_RF;
[comb_orig, comb_nans, n_comb_nans] = nansum_alt(comb_pre_sum);
comb_orig(comb_nans) = NaN;
% determine the expected standard deviation at each frequency:
std_comb_orig = sqrt(proc_expected_f_avgd_std^2 ./ nansum(SNR_resc_comb_RF.^2));
std_comb_orig(comb_nans) = NaN;

% normalize to expected st. dev. (note redefinition will not affect the normalized version):
norm_comb_spectrum = comb_orig ./ std_comb_orig;
R_comb_orig = 1 ./ std_comb_orig; % expected mean of 1-bin KSVZ axion in combined spectrum

if sum(isnan(comb_orig)) ~= n_comb_nans || sum(isnan(std_comb_orig)) ~= n_comb_nans
    warning('unexpected number of NaNs in combined spectrum');
end

% determine experimental mean and standard deviation:
comb_mean_experimental = nanmean(norm_comb_spectrum);
comb_std_experimental = nanstd(norm_comb_spectrum);

disp('spectra combined');

%% ============================================================================================== %%
%% PLOT COMBINED SPECTRUM
% plot the combined spectrum, unnormalized and with expected standard deviations and contibruing
% raw spectra counts:
if plot_comb_spectrum
    comb_abs_max = 1.02 * max(abs(comb_orig));
    axes_1 = [min_f_RF_GHz, max_f_RF_GHz, -comb_abs_max, comb_abs_max];
    
%     fig = figure('Position', [SD.default_left, SD.default_bot, SD.default_width * 1.15, ...
%         SD.default_height]);
    fig = figure('Position', [SD.default_left, SD.default_bot, SD.default_width * 2.0, ...
        1.1 * SD.default_height])
    hold on;
    
    % plot the number of contributing raw spectra:
    axes_2 = [axes_1(1), axes_1(2), -0.25, max(10,round(max(comb_cont_raw_spectra) * 1.05, -1))];
    
    yyaxis('right');
    p5 = plot(f_comb_GHz, comb_cont_raw_spectra, 'Color', SD.black, 'LineWidth', 0.75);
    p5.Color(4) = 0.15;
    
    axis(axes_2);
    ylabel('Contributing raw spectra');
    yyaxis('left');
    
    p1 = plot(f_comb_GHz, comb_orig, 'Color', SD.mygreen, 'LineWidth', 0.75);
    p2 = plot(f_comb_GHz, std_comb_orig, '-', 'Color', SD.myorange, 'LineWidth', 1);
    %    p3 = plot(f_comb_GHz, -std_comb_orig, '-', 'Color', SD.myorange, 'LineWidth', 1);
    lgd_ord = [p1, p2];
    lgd_labels = {'combined spectrum', '$\pm\sigma^\epsilon$'};
    axis(axes_1);
    xlabel('$\nu_\mathrm{RF}$ (GHz)');
    ylabel('Combined rescaled power excess');
    % title('Combined Spectrum');
    
    lgd_ord = [lgd_ord, p5];
    lgd_labels = [lgd_labels, '$n_\mathrm{raw}$'];
    mylegend(lgd_ord, lgd_labels, 'Location', 'SouthEast');
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Combined Spectrum.jpg'], 'jpeg');
        disp('combined spectrum figure saved successfully');
    end
end

%% ============================================================================================== %%
%% PLOT, HISTOGRAM NORMALIZED COMBINED SPECTRUM
% plot and histogram the normalized version of the combined spectrum. Typically, the combined
% spectrum has a flat normalziation for its center, though at least for DEADMAU5 that center is
% dwarfed in bandwidth by the wings, where the standard deviation grows larger. The normalzied
% combined spectrum should be normally distributed as long as the CLT holds and there is no axion
% (or it does not occupy many bins). Even if the CLT does not hold for the processed spectra, the
% enhancement in number of added spectra by a factor of the number of raw spectra (minus some IF
% cuts for some bins) should make it valid here:
if plot_and_hist_norm_comb_spectrum
    % main plot with histogram:
    % form corresponding lines to plot:
    sigma_experimental_line_lower = comb_mean_experimental - comb_std_experimental * ...
        ones(1, n_bins_comb);
    sigma_experimental_line_upper = comb_mean_experimental + comb_std_experimental * ...
        ones(1, n_bins_comb);
    
    % set axes for the combined  spectrum vs. frequency plot:
    max_y = 1.05 * max(abs(norm_comb_spectrum));
    axes_comb_spec = [min(f_comb_GHz), max(f_comb_GHz), -max_y, max_y];
    
    u_bound_edge = 1.00000001 * max(abs(norm_comb_spectrum)); % just above max abs power
    l_bound_edge = -u_bound_edge; % symmetric histogram
    bins_per_sigma = 30; % set a fixed number of bins per standard deviation
    % because the spectrum is normalized, its expected standard deviation is 1
    n_bins = round(bins_per_sigma * u_bound_edge);
    n_edges = n_bins + 1; % 1 more edge than bin
    edges = linspace(l_bound_edge, u_bound_edge, n_edges); % locations of edges
    edge_spacing = edges(2) - edges(1); % since edges are uniformly spaced, this works
    bin_centers = linspace(l_bound_edge + edge_spacing / 2, u_bound_edge - edge_spacing / 2, ...
        n_bins);
    counts = histcounts(norm_comb_spectrum, edges); % get counts in each bin
    counts_err = sqrt(counts);
    % calculate the expected Gaussian from the pure statistical expectation (assuming no axion)
    expected_gaussian = (edge_spacing * sum(~isnan(norm_comb_spectrum)) / sqrt(2 * pi)) * ...
        gaussmf(bin_centers, [1, 0]);
    axes_hist = [-max(counts) * 0.012, max([counts, expected_gaussian]) * 1.05, -max_y, max_y];
    
    fig = figure('Position', [SD.default_left, SD.default_bot, SD.default_width * 1.85, ...
        1.1 * SD.default_height]);
    pos_subp1 = [0.0930, 0.1430, 0.5520, 0.7539];
    subp1 = subplot('Position', pos_subp1);
    hold on;
    p1 = plot(f_comb_GHz, norm_comb_spectrum, 'LineWidth', 0.75, 'Color', SD.mygreen);
    p2 = plot(f_comb_GHz, sigma_experimental_line_lower, 'LineWidth', 1.5, 'Color', SD.mypink);
    p3 = plot(f_comb_GHz, sigma_experimental_line_upper, 'LineWidth', 1.5, 'Color', SD.mypink);
    p4 = plot(f_comb_GHz, ones(1, n_bins_comb), '--', 'LineWidth', 1.5, 'Color', SD.black);
    p5 = plot(f_comb_GHz, -ones(1, n_bins_comb), '--', 'LineWidth', 1.5, 'Color', SD.black);
    lgd_ord = [p1, p2, p4];
    lgd_labels = {'combined spectrum', '$\pm \sigma_\mathrm{meas}$', '$\pm 1$'};
    axis(axes_comb_spec);
    % title('Normalized Combined Spectrum');
    xlabel('$\nu_\mathrm{RF}$ (GHz)');
    ylabel('Normalized combined power excess');
    mylegend(lgd_ord, lgd_labels, 'Location', 'SouthEast');
    hold off;
    
    subp2 = subplot('Position', [pos_subp1(1) + pos_subp1(3),  pos_subp1(2), ...
        pos_subp1(3) * 0.55, pos_subp1(4)]);
    hold on;
    
    p7 = errorbar(counts, bin_centers,[],[],counts_err,counts_err, '.', 'Color', SD.mygreen);
    p8 = plot(expected_gaussian, bin_centers, 'Color', SD.black, 'LineWidth', 1);
    % title('Power Distribtuion');
    xlabel('Counts');
    mylegend([p7, p8], {'data', 'expected $\mathcal{N}$'}, 'Location', 'SouthEast');
    yticks([]); % no ticks on y-axis
    axis(axes_hist);
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Normalized Combined Spectrum.jpg'], 'jpeg');
        disp('normalized combined spectrum figure saved successfully');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % large, logarithmic histogram (uses same data as the histogram attached to the spectrum plot):
    max_x = u_bound_edge * 1.05;
    max_y_log_hist = 1.05 * max(log10([counts, expected_gaussian]));
    log_hist_axes = [-max_x, max_x, -0.1, max_y_log_hist];
    
    fig = figure;
    hold on;
    p1 = plot(bin_centers, log10(counts), '.', 'Color', SD.mygreen, 'MarkerSize', 22);
    p2 = plot(bin_centers, log10(expected_gaussian), 'Color', SD.black, 'LineWidth', 1.75);
    axis(log_hist_axes);
    % title('Combined Spectrum Power Distribution');
    xlabel('normalized power excess');
    ylabel('counts');
    if ceil(max_x) <= 10
        xticks(-ceil(max_x):ceil(max_x));
    end
    yticks(0:5);
    yticklabels({'$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'});
    mylegend([p1, p2],  {'data', 'expected $\mathcal{N}$'});
    text(-0.9, 1, ['$\mu$ = ', num2str(round(comb_mean_experimental, 4)), newline, ...
        '$\sigma$ = ', num2str(round(comb_std_experimental, 4))], 'Color', SD.mygreen);
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Normalized Combined Spectrum Histogram.jpg'], 'jpeg');
        disp('normalized combined spectrum histogram figure saved successfully');
    end
end

%% ============================================================================================== %%
%% COMBINED AUTOCORRELATION
% the correlation lengths at this point will not be affected by what lineshape we choose to apply
% for the grand spectrum, so they will give an indication of how much the SG filters affect the
% correlations in the combined data (presuming the SG's are the dominant contributor)
if plot_comb_autocor
    disp('calcualting combined spectrum autocorrelation...');
    corr_lag_extent_ind = round(1E6 * corr_lag_extent_MHz / raw_res_Hz);
    [lag_vec_norm_comb, autocor_norm_comb] = FUNC_autocorr_fcn_v2(norm_comb_spectrum, ...
        corr_lag_extent_ind);
    lag_vec_norm_comb_freq_kHz = lag_vec_norm_comb * raw_res_Hz*1E-3;
    disp([tab, 'combined spectrum autocorrelation calculated']);
    
    % calculate the moving mean
    k_win = 101; % window size (what would be 2W + 1, if W were the length in one direction)
    mo_avg_comb_autocor = movmean(autocor_norm_comb, k_win);
    
    expected_max_corr_length_kHz = 2 * max(SG_raw_filt_W_kHz, SG_norm_filt_W_kHz); % see ref. [10]
    
    ax = [0, max(lag_vec_norm_comb_freq_kHz), min(autocor_norm_comb), max(autocor_norm_comb)];
    
    fig = figure;
    hold on;
    p1 = plot(lag_vec_norm_comb_freq_kHz, autocor_norm_comb, 'Color', SD.mygreen, ...
        'LineWidth', 0.75);
    p2 = plot(lag_vec_norm_comb_freq_kHz, mo_avg_comb_autocor, 'Color', SD.black, ...
        'LineWidth', 1.2);
    p3 = plot(expected_max_corr_length_kHz * [1, 1], ax(3:4), '--');
    p4 = plot(ax(1:2), [0, 0]);
    axis(ax);
    % title('normalized combined spectrum autocorrelation');
    xlabel('$f_\mathrm{lag}$ (kHz)');
    ylabel('autocorrelation');
    mylegend([p2, p3], {[num2str(k_win), '-bin moving avg'], 'SG expected extent'});
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Combined Spectrum Autocorrelation.jpg'], 'jpeg');
        disp('combined spectrum autocorrelation figure saved successfully');
    end
end

%% ============================================================================================== %%
%% REDEFINE COMBINED SPECTRUM
% towards the first step in calculating the rebinned spectrum, ref. [2] actually (and somewhat
% annoyingly) redefines the combined spectrum and hence their standard deviations. I make explict
% note of these redefinitions rather than overwriting the same variable names, and since I am no
% longer rebinning, these will feed directly into the grand spectrum
disp(['redefining combined spectrum using Kg = ', num2str(Kg)]);
comb_redef = comb_orig * Kg;
std_comb_redef = std_comb_orig * Kg;

% ref. [2] uses D's with standard deviation R's where it once used deltas and sigmas. The D's are
% just the weighted excessres, where the ML weights are inverse-varainces
D_comb_RF = comb_redef ./ (std_comb_redef.^2);
% the R array is simply the inverse of the std array. It now (rescaled) represents the power
% deposited into a single bin by an axion that puts 1/Kg of its power there
R_comb_RF = 1 ./ std_comb_redef;
disp([tab, 'combined spectrum redefined']);

%% ============================================================================================== %%
%% AXION LINESHAPE
% the axion lineshape comes from ref. [2], with a full derivation/verification in ref. [13].
% Assumptions are made explicit there. The CDF is in ref. [15]

% write the pdf as an anonymous function - units of nu and nu_a do not matter as long as they are
% the same (e.g. GHz)
pdf_ax = @(nu, nu_a) (sqrt(6/pi) / (nu_a * r * EBeta2) * ...
    exp(-3 * (nu-nu_a) / (nu_a*EBeta2) - (3/2) * r^2) .* ...
    sinh(3 * r * sqrt(2 * (nu-nu_a) / (nu_a*EBeta2)))) .* ...
    heaviside(nu-nu_a); % the step function makes the PDF 0 below the axion rest frequency

% for the CDF, the inputs to the error functions must be real. Therefore I use absolute values for
% their nu-nu_a;s, with no harmul effects because the whole expression is multiplied by a step
% function
cdf_ax = @(nu, nu_a) (1/2) * (erf(sqrt(3/2)*r + sqrt(3*(abs(nu-nu_a))/(nu_a*EBeta2))) - ...
    erf(sqrt(3/2)*r - sqrt(3*(abs(nu-nu_a))/(nu_a*EBeta2))) - ...
    (sqrt(8/(3*pi))/r * exp(-(3*r^2)/2 - 3*(nu-nu_a)/(nu_a*EBeta2)) .* ...
    sinh(r*3*sqrt(2*(nu-nu_a)/(nu_a*EBeta2))))) .* heaviside(nu-nu_a);

disp('axion PDF and CDF functions defined');

%% ============================================================================================== %%
%% PLOT AXION PDF
if plot_axion_PDF
    axspec_discrete_spacing_kHz = 0.1; % just for plotting here
    axspec_dets_GHz = (ax_int_show_det_min_kHz:axspec_discrete_spacing_kHz:ax_int_kHz)*1E-6;
    n_axspec_dets = length(axspec_dets_GHz);
    axspec_freqs_GHz = ax_f_typ_GHz + axspec_dets_GHz;
    axspec_PDF_vals = pdf_ax(axspec_freqs_GHz, ax_f_typ_GHz);
    axspec_CDF_vals = cdf_ax(axspec_freqs_GHz, ax_f_typ_GHz); % use CDF formula
    % obtain the CDF by integration as a check:
    axspec_CDF_by_int_vals = NaN(1, n_axspec_dets);
    for i = 1:n_axspec_dets
        axspec_CDF_by_int_vals(i) = integral(@(nu)pdf_ax(nu, ax_f_typ_GHz), ...
            ax_f_typ_GHz, axspec_freqs_GHz(i));
    end
    
    % perform some mathematical checks
    ax_spec_int = integral(@(nu)pdf_ax(nu, ax_f_typ_GHz), ax_int_min_GHz, ...
        ax_int_max_GHz);
    ax_spec_trapz = trapz(axspec_freqs_GHz, axspec_PDF_vals);
    disp(['axion spectrum PDF integral = ', num2str(ax_spec_int, 4), ' (should be 1)']);
    disp(['axion spectrum PDF trapezoid sum = ', num2str(ax_spec_trapz, 4), ' (should be 1)']);
    min_x = min(axspec_freqs_GHz); max_x = max(axspec_freqs_GHz);
    min_y = 0; max_y = 1.02 * max(axspec_PDF_vals);
    
    sig_col = SD.myred;
    PDF_col = SD.myblue;
    CDF_col = SD.mygreen;
    
    fig = figure();
    hold on;
    
    % dots show samplings reflective of our experimental resolution
    p1 = plot(axspec_freqs_GHz, axspec_PDF_vals, '.', 'Color', PDF_col, 'MarkerSize', 16);
    p2 = plot(axspec_freqs_GHz, axspec_PDF_vals, '-', 'Color', PDF_col, 'LineWidth', 3);
    p3 = plot(ax_mean_GHz, min_y + 0.15*(max_y - min_y), '.', 'Color', sig_col, 'MarkerSize', 25);
    p4 = plot(ax_mean_GHz + [0, ax_lw_kHz*1E-6], (min_y + 0.15*(max_y - min_y)) * [1,1], '-', ...
        'Color', sig_col, 'LineWidth', 3 );
    
    xlim([min_x, max_x]);
    ylim([min_y, max_y]);
    % title(['$', num2str(ax_f_typ_GHz), '$ GHz axion, $' num2str(ax_int_kHz), '$ kHz int']);
    text(min_x + 0.6 * (max_x - min_x), min_y + 0.9 * (max_y - min_y), ...
        ['integral: ', num2str(ax_spec_int, 4)], 'Color', PDF_col);
    text(min_x + 0.6 * (max_x - min_x), min_y + 0.84 * (max_y - min_y), ...
        ['trapz: ', num2str(ax_spec_trapz, 4)], 'Color', PDF_col);
    text(ax_mean_GHz + ax_lw_kHz*1E-6, min_y + 0.15 * (max_y - min_y), ...
        '$\sigma_\mathrm{ax}$', 'Color', sig_col);
    text(ax_mean_GHz - 0.58 * ax_lw_kHz*1E-6, min_y + 0.15 * (max_y - min_y), ...
        '$\mu_\mathrm{ax}$', 'Color', sig_col);
    xlabel('$f$ (GHz)');
    ylabel('PDF');
    
    % plot CDF
    yyaxis('right');
    % plot the CDF obtained by integration and obtained directly
    p5 = plot(axspec_freqs_GHz, axspec_CDF_by_int_vals, '-', 'Color', CDF_col, 'LineWidth', 3);
    p6 = plot(axspec_freqs_GHz, axspec_CDF_vals, '.', 'MarkerSize', 16, 'Color', CDF_col);
    ylabel('CDF');
    yyaxis('left');
    plt = gca;
    plt.YAxis(1).Color = PDF_col;
    plt.YAxis(2).Color = CDF_col;
    
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Axion PDF.jpg'], 'jpeg');
        disp('axion PDF figure saved successfully');
    end
end

%% ============================================================================================== %%
%% COMBINED-GRAND WEIGHTS
% the grand spectrum is obtained by making a ML-weighted average of adjacent bins in the combined
% spectrum that takes into account the axion shape
disp('constructing combined-grand spectrum bin weights...');

% construct Lq from ref. [15], a modified version of that from ref. [2]. We will be testing many
% (correlated) hypotheses that the axion lies at specific frequencies, very closely spaced. Ben's
% approach was to test the hypothesis that it lay between 2 frequencies. This is a more
% straightforward idea (speaking strictly relatively)
% Lq says, up to a constant, what fraction of axions are in each contributing comb bin
Lq = NaN(grand_density, Kg);
disp([tab, 'constructing Lq (axion densitiy in each CS bin for each GS bin)...']);
% vector of possible misalignments of GS bins from CS bins (excepting left edge of GS) - negative
% implies that the axion is to the left of the CS bin center, 0 that it is at center, etc...
misalgn_delta_f_GHz = -raw_res_Hz*1E-9/2:grand_spacing_GHz:raw_res_Hz*1E-9/2;
misalgn_delta_f_GHz = misalgn_delta_f_GHz(1:end-1); % last entry is redundant with first
ax_f_misalgn_GHz = ax_f_typ_GHz + misalgn_delta_f_GHz;
for q = 1:Kg
    % frequency of nearest CS bin left edge - we imagine a CS bin centered right at the "typical"
    % axion. Then in the inner loop we will consider the CDF values for all allowed misalignments
    f_comb_left_GHz = ax_f_typ_GHz - raw_res_Hz*1E-9/2 + (q-1)*raw_res_Hz*1E-9;
    f_comb_right_GHz = f_comb_left_GHz + raw_res_Hz*1E-9;
    
    for j = 1:grand_density
        % if a frequency at which we are calculating the CDF is out of the integration window, reset
        % it to the edge of the window. This is a largely stylistic choice, but it enforces the
        % logic that we only look at the axion spectrum out to the edge of the integration window,
        % and no further, simplifying expectations and debugging
        f_eff_comb_left_GHz = min(ax_f_misalgn_GHz(j) + ax_int_kHz*1E-6, f_comb_left_GHz);
        f_eff_comb_right_GHz = min(ax_f_misalgn_GHz(j) + ax_int_kHz*1E-6, f_comb_right_GHz);
        
        % determine the value of the CDF at the start and end of the relevant combined bin:
        cdf_comb_left = cdf_ax(f_eff_comb_left_GHz, ax_f_misalgn_GHz(j));
        cdf_comb_right = cdf_ax(f_eff_comb_right_GHz, ax_f_misalgn_GHz(j));
        
        % Lq quantifies (up to Kg) the fraction of axion particles (hence power) in each combined
        % spectrum bin. Here, Lq id defined via 2 indices: one that tells which contibuting combined
        % bin (the 1st, 2nd, ..., Kg^th) that we are inquiring about. The second index is for how
        % far the axion we are considering is from the (above) 1st bin center (it will always be
        % within half a bin on either side). This index maps to the vector 'misalgn_delta_f_GHz',
        % where negative numbers mean an axion to the left of the first contributing combined bin's
        % center
        Lq(j,q) = Kg * (cdf_comb_right - cdf_comb_left);
    end
end
if sum(sum(isnan(Lq))) > 0
    warning('NaN(s) in Lq');
end
disp([tab, 'Lq calculated']);

%% ============================================================================================== %%
%% PLOT COMBINED-GRAND WEIGHTS
if plot_Lq
    min_x_q = 1 - 25;
    max_x_q = Kg + 25;
    min_y = 0;
    max_y = max(max(Lq/Kg)) * 1.1;
    
    col_ax2left = SD.myblue;
    col_ax2right = SD.myorange;
    comb_col = SD.mygreen;
    
    fig = figure;
    hold on;
    for q = 1:Kg
        plot(q*[1,1], [min_y, max_y], 'Color', SD.mygreen, 'LineWidth', 0.5);
        plot(q-1/2*[1,1], [min_y, max_y/20], 'Color', SD.black, 'LineWidth', 0.5);
        plot(q+1/2*[1,1], [min_y, max_y/20], 'Color', SD.black, 'LineWidth', 0.5);
    end
    for j = 1:grand_density
        if grand_density == 1
            col = col_ax2left;
        else
            col = (1-(j-1)/(grand_density-1))*col_ax2left + (j-1)/(grand_density-1)*col_ax2right;
        end
        
        plot(1/2 + (j-1)/grand_density * [1,1], [min_y, max_y], 'Color', col, 'LineWidth', 0.5);
        for q = 1:Kg
            plot([q-1/2,q+1/2], Lq(j,q)/Kg*[1,1], '-', 'LineWidth', 0.5, 'Color', col);
        end
        plot(1:Kg, Lq(j,:).'/Kg, '.', 'Color', col, 'MarkerSize', 12);
    end
    xlim([min_x_q, max_x_q]);
    ylim([min_y, max_y]);
    xlabel('$q$');
    ylabel('$L_q/K^g$');
    % title('$L_q$ and Misaligment');
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Lq and Misalignment.jpg'], 'jpeg');
        disp('Lq and misalignment figure saved successfully');
    end
end

%% ============================================================================================== %%
%% OBTAIN GRAND SPECTRUM
% first, map grand frequencies to combined ones, and then fill in the grand spectrum

disp([tab, 'mapping grand spectrum frequencies to combined spectrum frequencies...']);
% beacause of the lopsidedness of the grand spectrum, the first grand spectrum frequency can in
% principle be as far as we are going to integrate out on the axion  below the first combined
% spectrum bin (and then another half-bin, technically). The max is just that same
% technicality-half-bin above the last combined bin. Note: changing the start will necessitate
% changing how malign_ind is assigned, below
min_f_grand_GHz = min_f_RF_GHz - ax_int_kHz*1E-6 - raw_res_Hz*1E-9/2;
max_f_grand_GHz = max_f_RF_GHz + raw_res_Hz*1E-9/2;
f_grand_GHz = min_f_grand_GHz:grand_spacing_GHz:max_f_grand_GHz;
f_grand_hi_GHz = f_grand_GHz + ax_int_kHz*1E-6; % integrate grand spectrum bin to these end freqs
n_bins_grand = length(f_grand_GHz);

% we want the lowest combined freqeuency that contributes to each grand frequency. It is
% computationally quickest to calculate these all at once, but you need a for-loop to avoid storing
% a massive combined-by-grand matrix
g2c_lo_inds = NaN(1, n_bins_grand);
g2c_hi_inds = NaN(1, n_bins_grand);
disp([tab, 'mapping ', num2str(n_bins_grand), ' grand spectrum bins to ', num2str(n_bins_comb), ...
    ' combined spectrum start and stop bins...']);
g2c_lo_inds(1) = 1; g2c_hi_inds(1) = 1; % by construction, both mappings start at 1
eps = raw_res_Hz*1E-9*1E-5; % tiny number for numerical buffer
wb = waitbar(0, 'grand-comb mapping...');
for i = 2:n_bins_grand
    % a computationally faster way to compute the bin correspondances:
    % detuning to previously used combined frequency and to next (or final) one:
    del_to_prev_GHz = abs(f_grand_GHz(i) - f_comb_GHz(g2c_lo_inds(i-1)));
    del_to_next_GHz = abs(f_grand_GHz(i) - f_comb_GHz(min(g2c_lo_inds(i-1) + 1, n_bins_comb)));
    if del_to_prev_GHz < del_to_next_GHz - eps % if closer to previous one, use it again
        g2c_lo_inds(i) = g2c_lo_inds(i-1);
    else % use next one if not
        g2c_lo_inds(i) = min(g2c_lo_inds(i-1) + 1, n_bins_comb);
    end
    
    % similar for the high end:
    del_to_prev_GHz = abs(f_grand_hi_GHz(i) - f_comb_GHz(min(g2c_hi_inds(i-1), n_bins_comb)));
    del_to_next_GHz = abs(f_grand_hi_GHz(i) - f_comb_GHz(min(g2c_hi_inds(i-1) + 1, n_bins_comb)));
    if del_to_prev_GHz < del_to_next_GHz - eps
        g2c_hi_inds(i) = g2c_hi_inds(i-1);
    else
        g2c_hi_inds(i) = min(g2c_hi_inds(i-1) + 1, n_bins_comb);
    end
    
    if mod(i, 1E4) == 0
        waitbar(i/n_bins_grand, wb);
    end
end
close(wb);
if sum(isnan(g2c_lo_inds)) + sum(isnan(g2c_lo_inds)) > 0
    warning('NaN(s) in combined-grand mapping');
end
disp([tab, 'grand/combined mapping established']);


disp('constructing grand spectrum...');
R_grand_RF = NaN(1, n_bins_grand);
D_grand_RF = NaN(1, n_bins_grand);
wb = waitbar(0, 'filling grand spectrum...');
for l = 1:n_bins_grand
    % determine contribuing combined bins:
    start_ind = g2c_lo_inds(l);
    stop_ind = g2c_hi_inds(l);
    range = start_ind:stop_ind;
    Kg_delta = Kg - length(range);
    
    % specify the range of q's that get used (this is mostly just to handle edge cases, but seeing
    % the right sensitivity fall-off at the edges is a possible check
    q_start = 1;
    q_end = Kg;
    if Kg_delta > 0
        if start_ind == 1 % grand frequency is to left of left edge of left-most combined bin
            q_start = Kg - length(range) + 1; % start at a later q (only get right part of axion)
        elseif stop_ind == n_bins_comb % off the right edge
            q_end = length(range); % end at an earlier q (only get left part of axion)
        end
    end
    q_range = q_start:q_end;
    
    if length(range) > Kg
        error('range of contributing combined spectrum bins too large');
    end
    if length(q_range) ~= length(range)
        error('improper range for q in filling grand spectrum');
    end
    
    % to figure out which misalignment index to use, consider how the grand spectrum frequencies
    % were constructed: the first is one to the very left edge of a combined bin (misalignment index
    % 2), and the indices cyclically repeat with a frequency deptermined by the grand spectrum
    % density. Hence:
    malign_ind = mod(l-1, grand_density) + 1;
    
    % compute grand spectrum:
    D_grand_RF(l) = nansum_alt(D_comb_RF(range) .* Lq(malign_ind, q_range));
    R_grand_RF(l) = sqrt(nansum_alt((R_comb_RF(range) .* Lq(malign_ind, q_range)).^2));
    
    if mod(l, 1E4) == 0
        waitbar(l/n_bins_grand, wb);
    end
end
close(wb);
% this should have a Gaussian statistics, if purely noise
norm_grand_spectrum = D_grand_RF ./ R_grand_RF;
if sum(isnan(D_grand_RF)) + sum(isnan(R_grand_RF)) > 0
    warning('NaN(s) in grand spectrum');
end
if sum(isnan(norm_grand_spectrum)) > 0
    warning('NaN(s) in normalized grand spectrum');
end
disp([tab, 'grand spectrum constructed']);

% determine experimental mean and standard deviation:
grand_mean_experimental = nanmean(norm_grand_spectrum);
grand_std_experimental = nanstd(norm_grand_spectrum);

%% ============================================================================================== %%
%% GRAND STATISTICS
% how well does the grand spectrum match a Gaussian profile? A Gaussian is constrained by two
% parameters (DOF's), a mean and a standard deviation, and in chi2gof MATLAB estimates these from
% the data, and uses a test statistic (ref. [12]) to ask how likely (p-value) it is that we would
% observe that or more extreme given the data:
[~, p_val_all_data] = chi2gof(norm_grand_spectrum);

%% ============================================================================================== %%
%% PLOT, HISTOGRAM NORMALIZED GRAND SPECTRUM
if plot_and_hist_norm_grand_spectrum
    % main plot with histogram:
    % form corresponding lines to plot:
    sigma_experimental_line_lower = grand_mean_experimental - grand_std_experimental * ...
        ones(1, n_bins_grand);
    sigma_experimental_line_upper = grand_mean_experimental + grand_std_experimental * ...
        ones(1, n_bins_grand);
    
    % set axes for the rebinned spectrum vs. frequency plot:
    max_y = 1.05 * max(abs(norm_grand_spectrum));
    axes_grand_spec = [min_f_grand_GHz, max_f_grand_GHz, -max_y, max_y];
    
    u_bound_edge = 1.00000001 * max(abs(norm_grand_spectrum)); % just above max abs power
    l_bound_edge = -u_bound_edge; % symmetric histogram
    bins_per_sigma = 20; % set a fixed number of bins per standard deviation
    % because the spectrum is normalized, its expected standard deviation is 1
    n_bins = round(bins_per_sigma * u_bound_edge);
    n_edges = n_bins + 1; % 1 more edge than bin
    edges = linspace(l_bound_edge, u_bound_edge, n_edges); % locations of edges
    edge_spacing = edges(2) - edges(1); % since edges are uniformly spaced, this works
    bin_centers = linspace(l_bound_edge + edge_spacing / 2, u_bound_edge - edge_spacing / 2, ...
        n_bins);
    counts = histcounts(norm_grand_spectrum, edges); % get counts in each bin
    % calculate the naively (no correlations in the rebinned spectrum) expected Gaussian from the
    % pure statistical expectation (assuming no axion)
    expected_gaussian = (edge_spacing * sum(~isnan(norm_grand_spectrum)) / sqrt(2 * pi)) * ...
        gaussmf(bin_centers, [1, 0]);
    % define the estimated Gaussian to be that with the measured mean and standard deviation:
    estimated_gaussian = (edge_spacing * sum(~isnan(norm_grand_spectrum)) / (sqrt(2 * pi) * ...
        grand_std_experimental)) * gaussmf(bin_centers, [grand_std_experimental, ...
        grand_mean_experimental]);
    axes_hist = [-max(counts) * 0.012, max([counts, expected_gaussian]) * 1.05, -max_y, max_y];
    
    fig = figure('Position', [SD.default_left, SD.default_bot, SD.default_width * 1.85, ...
        1.1 * SD.default_height]);
    pos_subp1 = [0.0930, 0.1430, 0.5520, 0.7539];
    subp1 = subplot('Position', pos_subp1);
    hold on;
    p1 = plot(f_grand_GHz, norm_grand_spectrum, 'LineWidth', 1.5, 'Color', SD.mygreen);
    p2 = plot(f_grand_GHz, R_grand_RF, 'LineWidth', 1.5, 'Color', SD.saddlebrown);
    p3 = plot(f_grand_GHz, sigma_experimental_line_lower, 'LineWidth', 1.5, 'Color', SD.mypink);
    p4 = plot(f_grand_GHz, sigma_experimental_line_upper, 'LineWidth', 1.5, 'Color', SD.mypink);
    p5 = plot(f_grand_GHz, ones(1, n_bins_grand), '--', 'LineWidth', 1.5, 'Color', SD.black);
    p6 = plot(f_grand_GHz, -ones(1, n_bins_grand), '--', 'LineWidth', 1.5, 'Color', SD.black);
    lgd_ord = [p1, p2, p3, p5];
    lgd_labels = {'spectrum', ...
        '$\mu_\mathrm{KSVZ}^\mathrm{no\ corr}$' '$\pm 1 \sigma_\mathrm{meas}$', '$\pm 1$'};
    
    axis(axes_grand_spec);
    % title('Normalized Grand Spectrum');
    xlabel('$\nu_\mathrm{RF}$ (GHz)');
    ylabel('normalized power excess');
    mylegend(lgd_ord, lgd_labels, 'Location', 'SouthEast');
    hold off;
    
    subp2 = subplot('Position', [pos_subp1(1) + pos_subp1(3),  pos_subp1(2), ...
        pos_subp1(3) * 0.55, pos_subp1(4)]);
    hold on;
    p7 = plot(counts, bin_centers, '.', 'Color', SD.mygreen);
    p8 = plot(expected_gaussian, bin_centers, 'Color', SD.black, 'LineWidth', 1);
    p9 = plot(estimated_gaussian, bin_centers, 'Color', SD.mygreen, 'LineWidth', 1);
    % title('Power Distribtuion');
    xlabel('counts');
    mylegend([p7, p8, p9], {'data', 'na\"ive $\mathcal{N}$', ...
        'estimated $\mathcal{N}$'}, 'Location', 'SouthEast');
    yticks([]); % no ticks on y-axis
    axis(axes_hist);
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Normalized Grand Spectrum.jpg'], 'jpeg');
        disp('normalized grand spectrum figure saved successfully');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % large, logarithmic histogram (uses same data as the histogram attached to the spectrum plot):
    max_x = u_bound_edge * 1.05;
    max_y_log_hist = 1.05 * max(log10([counts, expected_gaussian]));
    log_hist_axes = [-max_x, max_x, -0.1, max_y_log_hist];
    
    fig = figure;
    hold on;
    p1 = plot(bin_centers, log10(counts), '.', 'Color', SD.mygreen, 'MarkerSize', 22);
    p2 = plot(bin_centers, log10(expected_gaussian), 'Color', SD.black, 'LineWidth', 1.75);
    p3 = plot(bin_centers, log10(estimated_gaussian), 'Color', SD.mygreen, 'LineWidth', 1.75);
    axis(log_hist_axes);
    % title('Grand Spectrum Power Distribution');
    xlabel('normalized power excess');
    ylabel('counts');
    if ceil(max_x) <= 10
        xticks(-ceil(max_x):ceil(max_x));
    end
    yticks(0:5);
    yticklabels({'$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'});
    mylegend([p1, p2, p3],  {'data', 'na\"ive $\mathcal{N}$', 'est. $\mathcal{N}$'});
    text(-0.9, 1, ['$\mu$ = ', num2str(round(grand_mean_experimental, 4)), newline, ...
        '$\sigma$ = ', num2str(round(grand_std_experimental, 4))], 'Color', SD.mygreen);
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Normalized Grand Spectrum Histogram.jpg'], 'jpeg');
        disp('normalized grand spectrum histogram figure saved successfully');
    end
end

%% ============================================================================================== %%
%% GRAND AUTOCORRELATION
if plot_grand_autocor
    disp('calcualting grand spectrum autocorrelation...');
    corr_lag_extent_ind = round(1E6 * corr_lag_extent_MHz / raw_res_Hz);
    [lag_vec_norm_grand, autocor_norm_grand] = ...
        FUNC_autocorr_fcn_v2(norm_grand_spectrum(~isnan(norm_grand_spectrum)), corr_lag_extent_ind);
    lag_vec_norm_grand_freq_kHz = lag_vec_norm_grand * raw_res_Hz*1E-3;
    disp([tab, 'grand spectrum autocorrelation calculated']);
    
    expected_max_corr_length_kHz = 2 * max(SG_raw_filt_W_kHz, SG_norm_filt_W_kHz); % see ref. [10]
    
    ax = [0, max(lag_vec_norm_grand_freq_kHz), min(autocor_norm_grand), max(autocor_norm_grand)];
    
    fig = figure;
    hold on;
    p1 = plot(lag_vec_norm_grand_freq_kHz, autocor_norm_grand, 'Color', SD.mygreen, 'LineWidth', 2.0);
    p2 = plot(expected_max_corr_length_kHz * [1, 1], ax(3:4), '--');
    p3 = plot(ax(1:2), [0, 0]);
    axis(ax);
    % title('normalized grand spectrum autocorrelation');
    xlabel('$f_\mathrm{lag}$ (kHz)');
    ylabel('autocorrelation');
    mylegend(p2, 'SG expected extent');
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Grand Spectrum Autocorrelation.jpg'], 'jpeg');
        disp('grand spectrum autocorrelation figure saved successfully');
    end
end

%% ============================================================================================== %%
%% CORRECT FOR GRAND SPECTRUM CORRELATIONS
% correct the grand spectrum for the variance-reduction that it experiences relative as a result of
% correlations in the combined spectrum

disp('renormalizing grand spectrum');
xi = grand_std_experimental / comb_std_experimental;
% ref. [2]'s delta^g_l/\tilde{simga}^g_l:
x = norm_grand_spectrum / xi; % has std = 1 by construction iff sigma_comb = 1
eta = eta_p / xi; % see ref. [2], Eq. (32)

% the renormalized grand spectrum has std. 1 and mean 0 (no axion), or mean given by the line below
% (KSVZ axion):
mu_KSVZ = eta * R_grand_RF;
% note we've switchted to a simple notation for the transition to analysis code. x is the excess,
% which should be a standard normal r.v. absent an axion. mu_KSVZ is the mean of a KSVZ axion whose
% rest mass is at the relevant grand spectrum bin

% determine experimental mean and standard deviation:
renorm_grand_mean_experimental = nanmean(x);
renorm_grand_std_experimental = nanstd(x); % FLAG - on 9/9/19, not coming out to 1

disp([tab, 'grand spectrum renormalized']);

%% ============================================================================================== %%
%% DERIVED QUANTITIES - ANALYSIS
% Ref. [21] underlies the math of the analysis, both frequentist and Bayesian.
desired_excl = 1 - desired_prior_red;
g_target_norm_vec = (g_target_min:g_target_step:g_target_max).';
n_gs = length(g_target_norm_vec);

% turn the cutoff percentages into indecies/frequencies. Err towards the edges, but keep in bounds
low_cutoff_freq_ind = max([floor(n_bins_grand * low_cutoff_freq_pct / 100), 1]);
hi_cutoff_freq_ind = min([ceil(n_bins_grand * hi_cutoff_freq_pct / 100), n_bins_grand]);
low_cutoff_freq_GHz = f_grand_GHz(low_cutoff_freq_ind);
hi_cutoff_freq_GHz = f_grand_GHz(hi_cutoff_freq_ind);
n_freqs_btw_cutoffs = hi_cutoff_freq_ind - low_cutoff_freq_ind + 1;

% FLAG - give more consideration to this value
n_adj_corr = grand_density * (Kg-1);

% frequentism:
Fn = desired_prior_red; % total single-bin false negative rate
% assuming 2 idential scans:
fn = 1 - sqrt(1 - Fn); % single-scan, single-bin false negative rate
xT = mu_target + sqrt(2)*erfinv(2*fn - 1); % threshold
fp = erfc(xT/sqrt(2)) / 2; % single-scan, single-bin false-positive rate
Fp = fp^2;  % total single-bin false-positive rate

% rescans: 
rescan_cluster_within_bins = round(rescan_cluster_within_Hz/(grand_spacing_GHz*1E9));

%% ============================================================================================== %%
%% FREQUENTISM
disp('frequentist testing for rescan candidates');
disp([tab, 'total two-scan false negative rate Fn = ', num2str(Fn,2)]);
disp([tab, 'single-scan false negative rate fn = ', num2str(fn, 2)]);
disp([tab, 'total two-scan false positive rate Fp = ', num2str(Fp, 2)]);
disp([tab, 'single-scan false positive rate fp = ', num2str(fp, 2)]);
% determine the coupling at which to test in the frequintist framework:
% for all of the analysis, couplings are normalized to KSVZ
g_test = sqrt(mu_target ./ mu_KSVZ);
% perform the tests:
clicks = x > xT;
click_inds = find(clicks);
n_clicks = sum(clicks); % clicks will likely come in correlated clusters
n_clicks_expected = fp * n_bins_grand;
disp([tab, 'frequentist testing complete: ', num2str(n_clicks), ' clicks occurred vs. ', ...
    num2str(n_clicks_expected, 2), ...
    ' expected. Clicks are likely to turn up in correlated clusters']);

%% ============================================================================================== %%
%% LIST CLICKS
% for each fully contiguous group of clicks, find the frequency and excess of the highest grand
% spectrum powers. These are then put in a form easily accessable to the user
n_clusters = 1; % init 0 
new_cluster = 1; % init true
cluster_start_inds = NaN(1, n_clicks); % largest we could possibly need
cluster_end_inds = NaN(1, n_clicks); % largest we could possibly need

i = 1; % init 1 
for i = 1:n_clicks
    c_ind = click_inds(i);
    
    if new_cluster
        cluster_start_inds(n_clusters) = c_ind; % add a start index
    end
    
    % what we care about is the distance from the start of the cluster, otherwise we could
    % artificially daisy-chain candidates into one
    if (i < n_clicks) && (click_inds(i+1) <= ...
            cluster_start_inds(n_clusters) + rescan_cluster_within_bins)
        new_cluster = 0; % next one part of same cluster
    else 
        cluster_end_inds(n_clusters) = c_ind; % add an end index
        if i < n_clicks
            n_clusters = n_clusters + 1;
            new_cluster = 1; % next one starts a new cluster
        end
    end
end
cluster_start_inds = cluster_start_inds(1:n_clusters);
cluster_end_inds = cluster_end_inds(1:n_clusters);

click_xmax = NaN(1, n_clusters);
click_xmax_f_GHz = NaN(1, n_clusters);
click_xmax_grand_inds = NaN(1, n_clusters);
click_xmax_muKSVZ = NaN(1, n_clusters); 
for i = 1:n_clusters
    cluster_range = cluster_start_inds(i):cluster_end_inds(i);
    f_cluster_GHz = f_grand_GHz(cluster_range); 
    x_cluster = x(cluster_range);
    muKSVZ_cluster = mu_KSVZ(cluster_range);
    [click_xmax(i), ind_max_cluster] = max(x_cluster); % find the highest excess within each cluster
    click_xmax_f_GHz(i) = f_cluster_GHz(ind_max_cluster);
    click_xmax_muKSVZ(i) = muKSVZ_cluster(ind_max_cluster);
    [~, click_xmax_grand_inds(i)] = min(abs(click_xmax_f_GHz(i) - f_grand_GHz));
end

%% ============================================================================================== %%
%% PLOT, HISTOGRAM RENORMALIZED GRAND SPECTRUM
if plot_and_hist_renorm_grand_spectrum
    % main plot with histogram:
    % form corresponding lines to plot:
    sigma_experimental_line_lower = renorm_grand_mean_experimental - ...
        renorm_grand_std_experimental * ones(1, n_bins_grand);
    sigma_experimental_line_upper = renorm_grand_mean_experimental + ...
        renorm_grand_std_experimental * ones(1, n_bins_grand);
    
    % set axes for the rebinned spectrum vs. frequency plot:
    max_y = 1.05 * max(abs(x));
    axes_grand_spec = [min_f_grand_GHz, max_f_grand_GHz, -max_y, max_y];
    
    u_bound_edge = 1.00000001 * max(abs(x)); % just above max abs power
    l_bound_edge = -u_bound_edge; % symmetric histogram
    bins_per_sigma = 20; % set a fixed number of bins per standard deviation
    % because the spectrum is normalized, its expected standard deviation is 1
    n_bins = round(bins_per_sigma * u_bound_edge);
    n_edges = n_bins + 1; % 1 more edge than bin
    edges = linspace(l_bound_edge, u_bound_edge, n_edges); % locations of edges
    edge_spacing = edges(2) - edges(1); % since edges are uniformly spaced, this works
    bin_centers = linspace(l_bound_edge + edge_spacing / 2, u_bound_edge - edge_spacing / 2, ...
        n_bins);
    counts = histcounts(x, edges); % get counts in each bin
    % calculate the naively (no correlations in the rebinned spectrum) expected Gaussian from the
    % pure statistical expectation (assuming no axion)
    expected_gaussian = (edge_spacing * sum(~isnan(x)) / sqrt(2 * pi)) * ...
        gaussmf(bin_centers, [1, 0]);
    % define the estimated Gaussian to be that with the measured mean and standard deviation:
    estimated_gaussian = (edge_spacing * sum(~isnan(x)) / (sqrt(2 * pi) * ...
        renorm_grand_std_experimental)) * gaussmf(bin_centers, [renorm_grand_std_experimental, ...
        renorm_grand_mean_experimental]);
    axes_hist = [-max(counts) * 0.012, max([counts, expected_gaussian]) * 1.05, -max_y, max_y];
    
    fig = figure('Position', [SD.default_left, SD.default_bot, SD.default_width * 1.85, ...
        1.1 * SD.default_height]);
    %pos_subp1 = [0.0930, 0.1430, 0.5520, 0.7539];
    %subp1 = subplot('Position', pos_subp1);
    hold on;
    p1 = plot(f_grand_GHz, x, 'LineWidth', 1.5, 'Color', SD.mygreen);
    p2 = plot(f_grand_GHz, sigma_experimental_line_lower, 'LineWidth', 1.5, 'Color', SD.mypink);
    p3 = plot(f_grand_GHz, sigma_experimental_line_upper, 'LineWidth', 1.5, 'Color', SD.mypink);
    p4 = plot(f_grand_GHz, ones(1, n_bins_grand), '--', 'LineWidth', 1.5, 'Color', SD.black);
    p5 = plot(f_grand_GHz, -ones(1, n_bins_grand), '--', 'LineWidth', 1.5, 'Color', SD.black);
    p6 = plot(f_grand_GHz, mu_KSVZ, 'Color', SD.saddlebrown);
    % p7 = plot(click_xmax_f_GHz, click_xmax, '.');
    lgd_ord = [p1, p2, p4, p6];
    lgd_labels = {'grand spectrum', '$\pm \sigma_\mathrm{meas}$', '$\pm 1$', '$\mu_\mathrm{KSVZ}$'};
    
    axis(axes_grand_spec);
    % title('Renormalized Grand Spectrum');
    xlabel('$\nu_\mathrm{RF}$ (GHz)');
    ylabel('Normalized power excess');
    mylegend(lgd_ord, lgd_labels, 'Location', 'SouthEast');
    hold off;
    
%     subp2 = subplot('Position', [pos_subp1(1) + pos_subp1(3),  pos_subp1(2), ...
%         pos_subp1(3) * 0.55, pos_subp1(4)]);
%     hold on;
%     p7 = plot(counts, bin_centers, '.', 'Color', SD.mygreen);
%     p8 = plot(expected_gaussian, bin_centers, 'Color', SD.black, 'LineWidth', 1);
%     p9 = plot(estimated_gaussian, bin_centers, 'Color', SD.mygreen, 'LineWidth', 1);
%     % title('Power Distribtuion');
%     xlabel('counts');
%     mylegend([p7, p8, p9], {'data', 'na\"ive $\mathcal{N}$', ...
%         'estimated $\mathcal{N}$'}, 'Location', 'SouthEast');
%     yticks([]); % no ticks on y-axis
%     axis(axes_hist);
%     hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Renormalized Grand Spectrum.jpg'], 'jpeg');
        disp('renormalized grand spectrum figure saved successfully');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % large, logarithmic histogram (uses same data as the histogram attached to the spectrum plot):
    max_x = u_bound_edge * 1.05;
    max_y_log_hist = 1.05 * max(log10([counts, expected_gaussian]));
    log_hist_axes = [-max_x, max_x, -0.1, max_y_log_hist];
    
    % include assumed-Poisson error bars:
    u_error = log10(1 + 1./sqrt(counts));
    l_error = log10(1 - 1./sqrt(counts));
    l_error(abs(l_error) == Inf) = 0; % dont show lower error bars that go down indefinitely
    fig = figure;
    hold on;
    p1 = errorbar(bin_centers, log10(counts), l_error, u_error, '.', 'Color', SD.mygreen, ...
        'MarkerSize', 22);
    p2 = plot(bin_centers, log10(expected_gaussian), 'Color', SD.black, 'LineWidth', 1.75);
    p3 = plot(bin_centers, log10(estimated_gaussian), 'Color', SD.mygreen, 'LineWidth', 1.75);
    axis(log_hist_axes);
    % title('Renormalized Grand Spectrum Power Distribution');
    xlabel('normalized power excess');
    ylabel('counts');
    if ceil(max_x) <= 10
        xticks(-ceil(max_x):ceil(max_x));
    end
    yticks(0:5);
    yticklabels({'$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'});
    mylegend([p1, p2, p3],  {'data', 'na\"ive $\mathcal{N}$', 'est. $\mathcal{N}$'});
    text(-0.9, 1, ['$\mu$ = ', num2str(round(renorm_grand_mean_experimental, 4)), newline, ...
        '$\sigma$ = ', num2str(round(renorm_grand_std_experimental, 4))], 'Color', SD.mygreen);
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'Renormalized Grand Spectrum Histogram.jpg'], 'jpeg');
        disp('renormalized grand spectrum histogram figure saved successfully');
    end
end

%% ============================================================================================== %%
%% UPDATE PRIORS & AGGREGATE
% update the priors as in ref. 21. Priors assumed unfirom over our (narrow) freqeuncy window
disp(['updating priors for ', num2str(n_gs), ' couplings and ', num2str(n_bins_grand), ...
    ' grand spectrum frequencies']);
mu = g_target_norm_vec.^2 * mu_KSVZ;
% define the function used to update priors (assumes sigma of the GS is 1):
u_func = @(x_meas, mu_ax) exp(-mu_ax.^2/2 + mu_ax .* x_meas);
u = u_func(x, mu); % apply the function to our data
u(isnan(u)) = 1; % set NaNs to unit update
disp([tab, 'priors updated']);

% calculate the corresponging aggreagte priors within the specified range:
disp('aggregating prior updates');

agg_u = mean(u(:, low_cutoff_freq_ind:hi_cutoff_freq_ind), 2);
disp([tab, 'prior updates aggregated']);

%% ============================================================================================== %%
%% ID CANDIATES
% idenfity the top several axion candidates within the cutoff range. Only one candidate per
% correlation-window centered frequency. Order from biggest candidates to smallest
cand_f_inds = NaN(1, n_cands);
cand_g_inds = NaN(1, n_cands);
cand_u = NaN(1, n_cands); % in the exact bin
cand_ext_u = NaN(1, n_cands); % extended over the correlation window

% make a reduced version of the total prior updates. This version will have candidates written over
% with 1's as candidates are identified
candID_u = u(:, low_cutoff_freq_ind:hi_cutoff_freq_ind);
for c = 1:n_cands
    [cand_u(c), cand_cutoff_f_ind] = max(max(candID_u));
    [~, cand_g_inds(c)] = max(max(candID_u.'));
    % the frequency index needs to be adjusted for the fact that it was in the structure that was
    % cut off at the ends:
    cand_f_inds(c) = cand_cutoff_f_ind + low_cutoff_freq_ind - 1;
    
    % take all couplings within the correlation window of the candidate and set their prior updates
    % to 1 (this seems more conservative than setting to 0, a nearby value, etc...):
    % handle edge cases by rounding to the first or last bin index
    lower_cand_corr_bound = max([1, cand_cutoff_f_ind - n_adj_corr]);
    upper_cand_corr_bound = min([n_freqs_btw_cutoffs, cand_cutoff_f_ind + n_adj_corr]);
    
    % also use these bounds to set the extended prior updates. Doing this with the candidate-ID'ing
    % structure (before resetting the relevant updates to 1) means we will always count overlaps
    % towards the bigger candidate, and not the smaller
    cand_ext_u(c) = sum(candID_u(cand_g_inds(c), ...
        lower_cand_corr_bound:upper_cand_corr_bound));
    
    candID_u(:, lower_cand_corr_bound:upper_cand_corr_bound) = 1;
end

%% ============================================================================================== %%
%% FREQUENTIST AGGREGATION
% calculate the single-scan false negative rate:
f_neg_spec = 1/2 * (erf((-mu_target * (g_target_norm_vec.^2 ./ g_test.^2) + xT) / (sqrt(2))) + 1);
% total false negative rate for an n-scan threshold-based protocol:
f_neg_spec = 1 - (1 - f_neg_spec).^n_thresh_scans;

% replace all NaN's with unit false negative rates. These inevitably bring up the aggregate
% exclusion to stronger values, just as in the case of the aggregate prior updates. Good for an
% apples-to-apples comparison:
% FLAG - make sure this is behaving properly
f_neg_spec(isnan(f_neg_spec)) = 1;

% define an aggregate exlcusion over the specified range (see ref. [2]) that comports with the
% aggregate prior reduction:
agg_freq_excl = 1 - mean(f_neg_spec(:, low_cutoff_freq_ind:hi_cutoff_freq_ind), 2);

%% ============================================================================================== %%
%% SUB-AGGREGATE
% for the color exclusion plot, we want to cut down on the number of points displayed in the plots
% (see ref. [4]). I have removed the old "smoothing" that I had, as that entails averaging in a bad
% space
n_red_edges = n_red_freqs + 1; % have one more edge to the reduced frequency bins than actual bins
eps = 1E-12; % much smaller than spacing between freq bins
red_edges_GHz = linspace(min_f_grand_GHz - eps, max_f_grand_GHz + eps, n_red_edges);
red_freq_spacing_GHz = red_edges_GHz(2) - red_edges_GHz(1);
red_freq_min_GHz = mean(red_edges_GHz(1:2));
red_freq_max_GHz = mean(red_edges_GHz(end - 1:end));
red_freqs_GHz = linspace(red_freq_min_GHz, red_freq_max_GHz, n_red_freqs);
% now make a vector the same length as the main frequency vector that maps each entry to the index
% of the reduced frequency vector. As discussed in ref. [4], simply map to the closest bin:
main2red_freq_map = NaN(1, n_bins_grand);
for i = 1:n_bins_grand
    [~, main2red_freq_map(i)] = min(abs(red_freqs_GHz - f_grand_GHz(i)));
end

% now we have to sub-aggreagte everything that appears on the color plot. This includes the Bayesian
% prior reduction line, the frequentist exclusion line, and the color map itself. Don't bother
% sub-aggregating things that will not be plotted

% sub-aggregate the color map:
sub_agg_u = NaN(n_gs, n_red_freqs);
for i = 1:n_red_freqs
    sub_agg_u(:, i) = mean(u(:, main2red_freq_map == i), 2);
end

% with the sub-aggreagted map in hand, now ask: what is the smallest coupling, if any, at which we
% get at least (to be conservative) the desired prior reduction? This will sub-aggregate the
% Bayesian "exclusion" line:
sub_agg_Bayes_excl_line = nan(1, n_red_freqs);
for j = 1:n_red_freqs
    i = 0;
    while i < n_gs && sub_agg_u(end - i, j) <= desired_prior_red
        sub_agg_Bayes_excl_line(j) = g_target_norm_vec(end - i);
        i = i + 1; % increment coupling index
    end
end

% lastly, sub-aggregate the frequensist exlcusion line. Do this by making a full map, and then,
% analagously to the Bayesian case, extracting the desired line from it intelligently:
sub_agg_freq_excl = nan(n_gs, n_red_freqs);
for i = 1:n_red_freqs
    sub_agg_freq_excl(:, i) = 1 - mean(f_neg_spec(:, main2red_freq_map == i), 2);
end
% find the smallest coupling, if any, at each frequnecy, at which we get at least the desired
% exclusion (analagous to above with an inequality reversal)
sub_agg_freq_excl_line = nan(1, n_red_freqs);
for j = 1:n_red_freqs
    i = 0;
    while i < n_gs && sub_agg_freq_excl(end - i, j) >= desired_excl
        sub_agg_freq_excl_line(j) = g_target_norm_vec(end - i);
        i = i + 1; % increment coupling index
    end
end

% mark the rescans along this line (put the frequencies at the center of sub-aggregate bins, else it
% looks weird):
sub_agg_rescan_freqs = red_freqs_GHz(main2red_freq_map(click_inds));
sug_agg_rescan_vals = sub_agg_freq_excl_line(main2red_freq_map(click_inds));

%% ============================================================================================== %%
%% FLOOR + CEILING
% put data in a nice form for plotting. Prior updates often go down to absurdly low levels (a 10^-30
% update is somewhat silly on physical grounds - we would be prepared to believe some glitch in our
% comptuer before betting at those odds!), so I will round a version of the total prior updates to
% reflect this. Do all this to the sub-aggregated data, which is what will be being plotted:
sub_agg_u_bound = max(sub_agg_u, u_flr);
sub_agg_u_bound = min(sub_agg_u_bound, u_ceil);
log_sub_agg_tot_u_bound = log10(sub_agg_u_bound);

% also do for the normal version of the data
tot_u_bound = max(u, u_flr);
tot_u_bound = min(tot_u_bound, u_ceil);
log_tot_u_bound = log10(tot_u_bound);

%% ============================================================================================== %%
%% CALCULATE SPEED-UP
% assuming a conclusion of exclusion, an important figure of merit of the BPM method is the time it
% takes to exclude. Since exclusion, or prior reduction, to a given certainty is not uniform in
% frequency, a single number must be computed. The aggregation does the bulk of the work here, but
% the last step is interpolating between the values to get the coupling at which we have a given
% amount of prior reduction (or the equivalent exclusion). To get the speed-up, scan time goes, all
% else equal, as coupling^(-4)
g_at_des_red = interp1(agg_u(agg_u > 0), g_target_norm_vec(agg_u > 0), desired_prior_red);
eps = 1E-6;
g_at_des_excl = interp1(agg_freq_excl(agg_freq_excl > eps & agg_freq_excl < 1-eps), ...
    g_target_norm_vec(agg_freq_excl > eps & agg_freq_excl < 1-eps), desired_excl);

% taking the exclusion coupling as the baseline to compare to, these are the speed-ups (values below
% 1 are slow-downs) of the other aggregate quantities at the desired exclusion/prior reduction
speed_up = (g_at_des_excl / g_at_des_red)^4; % figure of merit for BPM

% report the speed-ups
disp([num2str(100 * desired_excl), '% aggregate frequentist exclusion occurs at ', ...
    num2str(round(g_at_des_excl, 2)), ' KSVZ']);
disp(['equivalent prior reduction to ', num2str(100 * desired_prior_red), '% occurs at: ']);
disp(['    BPM: ', num2str(round(g_at_des_red, 2)), ' KSVZ (', num2str(round(speed_up, 2)), ...
    'X speed-up)']);

%% ============================================================================================== %%
%% SAVE EXCLUSION
% save outputs having to do with the exclusion to pass along to the rescan processing code
scan1excl.norm_cand_mks_szs_to_priors = norm_cand_mks_szs_to_priors;
scan1excl.min_f_RF_GHz = min_f_RF_GHz;
scan1excl.max_f_RF_GHz = max_f_RF_GHz;
scan1excl.n_red_freqs = n_red_freqs;
scan1excl.g_target_min = g_target_min;
scan1excl.g_target_max = g_target_max;
scan1excl.g_target_step = g_target_step;
scan1excl.g_target_show_top = g_target_show_top;
scan1excl.u_flr = u_flr;
scan1excl.u_ceil = u_ceil;
scan1excl.desired_prior_red = desired_prior_red;
scan1excl.low_cutoff_freq_pct = low_cutoff_freq_pct;
scan1excl.hi_cutoff_freq_pct = hi_cutoff_freq_pct;
scan1excl.n_cands = n_cands;
scan1excl.mu_target = mu_target;
scan1excl.n_thresh_scans = n_thresh_scans;
scan1excl.rescan_cluster_within_Hz = rescan_cluster_within_Hz;
scan1excl.click_xmax = click_xmax;
scan1excl.click_xmax_f_GHz = click_xmax_f_GHz;
scan1excl.click_xmax_grand_inds = click_xmax_grand_inds;
scan1excl.click_xmax_muKSVZ = click_xmax_muKSVZ;
scan1excl.mu = mu;
scan1excl.u = u;
scan1excl.agg_u = agg_u;
scan1excl.cand_f_inds = cand_f_inds;
scan1excl.cand_g_inds = cand_g_inds;
scan1excl.cand_u = cand_u;
scan1excl.cand_ext_u = cand_ext_u;
scan1excl.sub_agg_Bayes_excl_line = sub_agg_Bayes_excl_line;
scan1excl.mu_KSVZ = mu_KSVZ; 
scan1excl.f_neg_spec = f_neg_spec;
scan1excl.x = x;

if save_scan1excl_to_wkspc % for use in rescan processing
    disp('saving scan 1 exclusion data to workspace...');
    % create directory if necessary:
    if ~exist(save_to_folder, 'dir')
        mkdir(save_to_folder);
    end
    save([save_to_folder, f, save_to_scan1excl], varname(scan1excl), '-v7.3');
    disp([tab, 'scan 1 exlcusion data saved to workspace']);
end
    
%% ============================================================================================== %%
%% PLOT EXCLUSION
if plot_excl
    % aggregate prior versus axion-photon coupling:
    axes_1 = [g_target_min, g_target_show_top, 0, 1.01 * max(agg_u)];
    axes_2 = [g_target_min, g_target_show_top, (1 - 1.01 * max(agg_u)), 1];
    
    fig_w = 1400;
    fig_h = 600;
    new_fs = adj_f_sz_v1(fig_w / 2); % / 2 assumes that this plot will run across a whole page
    
    fig = figure('Position', [402, 84, fig_w, fig_h]); % left, bottom, width, height
    
    left_inc = -0.06;
    bot_inc = 0.065;
    width_mult = 1.1;
    height_mult = 0.9;
    
    % modify the things being plotted on the p-color plot to reflect the fact that MATLAB does not
    % show the right-most bin on its own (this is an incredibly stupid bug and I wish they would fix
    % it). This is technically also a problem for the top bin, but if we space couplings closely, as
    % we are free to do and will do for publication, the problem vanishes.
    log_sub_agg_u_bounded_padjust = [log_sub_agg_tot_u_bound, NaN(n_gs, 1)];
    % technically all the frequencies plotted correspond to left bin edges. The labels that are
    % affixed to them correspond (before being modified as below) to the bin centers. The left edges
    % are at half-a-bin lower frequencies
    red_freqs_padjust_GHz = [red_freqs_GHz, red_freq_max_GHz + red_freq_spacing_GHz] - ...
        red_freq_spacing_GHz / 2;
    
    % color exclusion plot:
    p_axes = [min_f_grand_GHz, max_f_grand_GHz, 0, g_target_show_top];
    % specify lines at the aggregation-cutoff frequencies:
    low_cutoff_line_x = [low_cutoff_freq_GHz, low_cutoff_freq_GHz];
    hi_cutoff_line_x = [hi_cutoff_freq_GHz, hi_cutoff_freq_GHz];
    low_hi_cutoff_line_y = [g_target_min, g_target_show_top];
    % make the color exlcusion plot:
    subp2 = subplot(1, 2, 1);
    hold on;
    
    pos_subp2 = get(subp2, 'Position'); % left, bottom, width, height
    set(subp2, 'Position', (pos_subp2 + [left_inc, bot_inc, 0, 0]) .* [1, 1, width_mult, ...
        height_mult]);
    h = pcolor(red_freqs_padjust_GHz, g_target_norm_vec, log_sub_agg_u_bounded_padjust);
    axis(p_axes);
    
    
    c = mycolorbar('Ticks', [log10(u_flr), log10(0.2), 0, log10(5), ...
        log10(u_ceil)], 'TickLabels', {'$\leq 0.05$', '0.2', '1', '5', '$\geq 20$'});
    set(c.Title, 'string', '$U_s$', 'interpreter', 'latex');
    
    label_text_frac_x = 0.035;
    label_text_x = p_axes(1) + label_text_frac_x * (p_axes(2) - p_axes(1));
    label_text_frac_y = 0.925;
    
    
    p1 = plot(red_freqs_GHz, sub_agg_Bayes_excl_line, 'Color', SD.myblue, 'LineWidth', 2);
    p2 = plot(red_freqs_GHz, sub_agg_freq_excl_line, 'Color', SD.myred, 'LineWidth', 2);
    % draw in the aggregation cutoff lines
    color_bounds = SD.black;
    p3 = plot(low_cutoff_line_x, low_hi_cutoff_line_y, '--', 'Color', color_bounds, 'LineWidth', 1.0);
    p4 = plot(hi_cutoff_line_x, low_hi_cutoff_line_y, '--', 'Color', color_bounds, 'LineWidth', 1.0);
    
    
    
    yticks(round(g_target_min):1:round(g_target_show_top));
    % mark the top candidates with markers that scale with the prior-update size of the candidate:
    pos_fig_pxl = get(fig, 'Position'); % left, bottom, width, height
    width_fig_pxl = pos_fig_pxl(3);
    pos_subp2_frac = get(subp2, 'Position'); % left, bottom, width, height
    width_subp2_pxl = width_fig_pxl * pos_subp2(3);
    % the width of the markers relative to the plot will represent the fraction of our priors or
    % posteriors take up at the coupling where their prior update maximizes:
    cand_linewidths = 1.5 * ones(1, n_cands);
    if norm_cand_mks_szs_to_priors % normalize to priors
        cand_mkr_sz = width_subp2_pxl * cand_ext_u / n_bins_grand;
    else % normalize to posteriors
        cand_mkr_sz = width_subp2_pxl * cand_ext_u ./ (agg_u(cand_g_inds).' * n_bins_grand);
    end
    for c = 1:n_cands
        cand_marker_col = SD.neoncarrot;
        
        marker_center = red_freqs_GHz(main2red_freq_map(cand_f_inds(c)));
        marker_extent = (max_f_grand_GHz - min_f_grand_GHz) * cand_mkr_sz(c) / width_subp2_pxl;
        marker_left = marker_center - marker_extent / 2;
        marker_right = marker_center + marker_extent / 2;
        
        p5c = plot(red_freqs_GHz(main2red_freq_map(cand_f_inds(c))), ...
            g_target_norm_vec(cand_g_inds(c)), '.', 'Color', SD.black, 'MarkerSize', 36, ...
            'LineWidth', cand_linewidths(c));
        p5c = plot(red_freqs_GHz(main2red_freq_map(cand_f_inds(c))), ...
            g_target_norm_vec(cand_g_inds(c)), '.', 'Color', cand_marker_col, 'MarkerSize', 30, ...
            'LineWidth', cand_linewidths(c));
    end
    
    
    if mark_rescan_freqs
        rescan_marker_col = SD.black;
        p6b = plot(sub_agg_rescan_freqs, sug_agg_rescan_vals, '.', 'Color', rescan_marker_col, ...
            'MarkerSize', 1.0 * SD.main_marker_sz);
    end
    
    
    % solidify a plot boundary:
    plot([min_f_grand_GHz, max_f_grand_GHz], g_target_show_top * [1, 1], 'LineWidth', 4);
    
    xlabel('$\nu\ \textrm{(GHz)}$');
    ylabel('$|g_\gamma / g_\gamma^\mathrm{KSVZ}|$');
    % ylabel('$|g_\gamma / g_\gamma^\mathrm{KSVZ}|$');
    % set(subp2, 'LineWidth', 4);
    
    plt_main = gca;
    if show_LU_inset
        % parameter plot boundaries:
        ins_x_min = min_f_grand_GHz;
        ins_x_max = max_f_grand_GHz;
        ins_y_min = min([(cand_ext_u / n_bins_grand) * 0.75, 0.005]);
        ins_y_max = max(cand_ext_u ./ (agg_u(cand_g_inds).' * n_bins_grand)) * 1.5;
        ins_axes = [ins_x_min, ins_x_max, log10(ins_y_min), log10(ins_y_max)];
                
        ins_pos = plt_main.Position; % left, bottom, width, height
        start_b_g = 0; % starting bottom of the inset plot - in units of axion_photon coupling
        height_b_g = g_target_min;
        frac_w = 1; % fractional width of the inset plot
        frac_h = height_b_g / g_target_show_top; % fractional height
        ins_pos(1) = ins_pos(1);
        ins_pos(2) = ins_pos(2);
        ins_pos(3) = frac_w * ins_pos(3);
        ins_pos(4) = frac_h * ins_pos(4);
        
        axes('Position', ins_pos);
        hold on;
        
        axis(ins_axes);
        xticks([]);
        yticks([-2, -1]);
        yticklabels({'$1\%$', '10\%'});
        
        for c = 1:n_cands
            bar_lw = 5.5;
            bar_center = red_freqs_GHz(main2red_freq_map(cand_f_inds(c)));
            bar_space_frac = 0.001;
            bar_left = bar_center - (ins_x_max - ins_x_min) * bar_lw * bar_space_frac;
            bar_right = bar_center + (ins_x_max - ins_x_min) * bar_lw * bar_space_frac;
            bar_top_prior = cand_ext_u(c) / n_bins_grand;
            log_bar_top_prior = log10(bar_top_prior);
            bar_side_lw = 0.75;
            
            plot(bar_center * [1, 1], [log10(ins_y_min), log_bar_top_prior], '-', ...
                'LineWidth', bar_lw, 'Color', cand_marker_col);
            plot(bar_left * [1, 1], [log10(ins_y_min), log_bar_top_prior], '-', ...
                'LineWidth', bar_side_lw, 'Color', SD.black);
            plot(bar_right * [1, 1], [log10(ins_y_min), log_bar_top_prior], '-', ...
                'LineWidth', bar_side_lw, 'Color', SD.black);
            plot([bar_left, bar_right], log_bar_top_prior * [1, 1], '-', ...
                'LineWidth', bar_side_lw, 'Color', SD.black);
        end
        plot([ins_x_min, ins_x_max], log10(ins_y_min) * [1, 1], 'LineWidth', 4);
        plot([ins_x_min, ins_x_max], log10(ins_y_max) * [1, 1], 'LineWidth', 4);
        
        plt = gca;
        plt.YAxis(1).Color = cand_marker_col;
        
        t_height = log10(0.05);
        fs = new_fs;
        text(label_text_x, t_height, '$P_\mathrm{LU}^{(1)}/\mathcal{P}_a$', ...
            'Color', cand_marker_col, 'FontSize', fs);
        hold off;
    end
    
    subp1 = subplot(1, 2, 2);
    pos_subp1 = get(subp1, 'Position'); % left, bottom, width, height
    set(subp1, 'Position', (pos_subp1 + [0.4 * left_inc, bot_inc, 0, 0]) .* [1, 1, width_mult, ...
        height_mult]);
    hold on;
    axis(axes_1);
    
    label_text_x = axes_1(1) + label_text_frac_x * (axes_1(2) - axes_1(1));
    label_text_y = axes_1(3) + label_text_frac_y * (axes_1(4) - axes_1(3));
    
    p2 = plot(g_target_norm_vec, agg_u, 'Color', SD.myblue, 'LineWidth', 3);
    p3 = plot(g_target_norm_vec, desired_prior_red * ones(size(g_target_norm_vec)), '--', ...
        'Color', SD.mygray, 'LineWidth', 3);
    
    %     title(['Aggregate Updates for ', num2str(low_cutoff_freq_GHz), '-', ...
    %         num2str(hi_cutoff_freq_GHz), ' GHz']);
    xlabel('$|g_\gamma / g_\gamma^\mathrm{KSVZ}|$');
    ylabel('$\mathcal{U}$');
    
    % make ticks along the desired prior reduction line:
    y_tick_extent = 0.03;
    y_tick_min = desired_prior_red - y_tick_extent / 2;
    y_tick_max = desired_prior_red + y_tick_extent / 2;
    taus_to_plot = [64, 32, 16, 8, 4, 2, 1, 0.5];
    n_ticks = length(taus_to_plot);
    x_tick_locs = g_at_des_excl * taus_to_plot .^ (-1/4);
    x_tick_subt = 0.05 * ones(size(x_tick_locs));
    x_tick_subt(1) = 0.13;
    x_tick_subt(2) = 0.09;
    x_tick_subt(3) = 0.06;
    % call tau the time it would take for a scan. Normalize tau = 1 to the case of 90% frequentist
    % exclusion
    % plot in the ticks, including textual labels:
    for i = 1:n_ticks
        tick_col = SD.mygray;
        plot([x_tick_locs(i), x_tick_locs(i)], [y_tick_min, y_tick_max], 'Color', tick_col, ...
            'LineWidth', 3);
        text(x_tick_locs(i) - x_tick_subt(i), y_tick_max + 1.2 * y_tick_extent, ...
            num2str(taus_to_plot(i), 2), 'Color', tick_col);
    end
    text(2.37, 0.06, '$\tau_\mathrm{rel}$', 'Color', tick_col);
    
    % plot the aggregate frequentist exclusion:
    yyaxis('right');
    p6 = plot(g_target_norm_vec,  agg_freq_excl, 'Color', SD.myred, 'LineWidth', 3);
    axis ij; % reverses the y-axis
    axis(axes_2);
    ylabel('$\mathcal{E}$');
    
    yyaxis('left');
    plt = gca;
    plt.YAxis(1).Color = SD.myblue; % change color of LHS y-axis
    plt.YAxis(2).Color = SD.myred; % change color of RHS y-axis
    
    % inset emphasizing speedup:
    if show_speedup_inset
        % parameter plot boundaries:
        ins_x_min = 1.35;
        ins_x_max = 1.65;
        ins_y_min = 0.06;
        ins_y_max = 0.25;
        ins_axes_1 = [ins_x_min, ins_x_max, ins_y_min, ins_y_max];
        ins_axes_2 = [ins_axes_1(1), ins_axes_1(2), 1 - ins_axes_1(4), 1 - ins_axes_1(3)];
        
        % draw a box on the main plot
        plot([ins_x_min, ins_x_max, ins_x_max, ins_x_min, ins_x_min], ...
            [ins_y_min, ins_y_min, ins_y_max, ins_y_max, ins_y_min], '--', 'LineWidth', 2);
        
        % define axes of the inset:
        ins_pos = plt.Position; % left, bottom, width, height
        frac_w = 0.4; % fractional width of the inset plot
        frac_h = 0.4; % fractional height of the inset plot
        ins_pos(1) = ins_pos(1) + 0.56 * ins_pos(3);
        ins_pos(2) = ins_pos(2) + 0.555 * ins_pos(4);
        ins_pos(3) = frac_w * ins_pos(3);
        ins_pos(4) = frac_h * ins_pos(4);
        
        axes('Position', ins_pos);
        hold on;
        
        axis(ins_axes_1);
        
        x_range = axes_1(2) - axes_1(1);
        ins_x_range = ins_axes_1(2) - ins_axes_1(1);
        frac_ins_range = ins_x_range / x_range;
        lw_expand_frac = 1.25;
        
        ins_p2 = plot(g_target_norm_vec, agg_u, 'Color', SD.myblue, 'LineWidth', 3*lw_expand_frac);
        ins_p3 = plot(g_target_norm_vec, desired_prior_red * ones(size(g_target_norm_vec)), ...
            '--', 'Color', SD.mygray, 'LineWidth', 3);
        
        % plot ticks within the inset at the 3 crossover points of interest
        ins_y_tick_extent = y_tick_extent * frac_h;
        ins_y_tick_min = desired_prior_red - ins_y_tick_extent / 2;
        ins_y_tick_max = desired_prior_red + ins_y_tick_extent / 2;
        
        ins_x_tick_locs = [g_at_des_red, g_at_des_excl];
        ins_taus_to_plot = (ins_x_tick_locs / g_at_des_excl) .^ -4;
        
        % plot the 2 ticks and labels manually:
        fs = 0.85 * new_fs;
        % 1. total aggregate prior reduction:
        plot([ins_x_tick_locs(1), ins_x_tick_locs(1)], [ins_y_tick_min, ins_y_tick_max], ...
            'Color', SD.myblue, 'LineWidth', 3);
        text(ins_x_tick_locs(1) - 0.14 * frac_h, ins_y_tick_min - 3.4 * ins_y_tick_extent / 3, ...
            num2str(ins_taus_to_plot(1), 3), 'Color', SD.myblue, 'FontSize', fs);
        % 2. aggregate exclusion:
        plot([ins_x_tick_locs(2), ins_x_tick_locs(2)], [ins_y_tick_min, ins_y_tick_max], ...
            'Color', SD.myred, 'LineWidth', 3);
        text(ins_x_tick_locs(2) - 0.03 * frac_h, ins_y_tick_max + 2.4 * ins_y_tick_extent / 3, ...
            num2str(ins_taus_to_plot(2), 3), 'Color', SD.myred, 'FontSize', fs);
        
        xticks([]);
        yticks([]);
        
        yyaxis('right');
        ins_p6 = plot(g_target_norm_vec,  agg_freq_excl, 'Color', SD.myred, ...
            'LineWidth', 3*lw_expand_frac);
        axis ij; % reverses the y-axis
        axis(ins_axes_2);
        yticks([]);
        hold off;
    end
    hold off;
    
    if save_figs
        saveas(fig, [save_to_folder, f, 'BPM Exclusion.jpg'], 'jpeg');
        disp('BPM aggregate exclusion and map figure saved successfully');
    end
end

%% ============================================================================================== %%
%% TIMING RESULTS
if print_timing
    dt_stop = datetime; % record stop time
    t_run_s = etime(datevec(dt_stop), datevec(dt_start)); % run time
    disp(['run time: ', num2str(round(t_run_s / 60, 2)), ' minutes']);
end

%% #################################################################################################
%% ######################################## END OF PROGRAM #########################################
%% #################################################################################################