%% DESCRIPTION
% Reads in data in the fashion of HAYSTACp2_process_v2, and uses that data to construct a set
% of workspaces containing only the necessary information for HAYSTAC data processing for a raw
% spectra. Each raw spectra is self-contained
%
%% HISTORY
% - v1 created by Dan Palken on 8/28/19. Goes with HAYSTACp2_process_v2
% - v2 created by DP on 8/30/19
% - v3 created by DP on 8/30/19. Takes in different frequency vectors for the different ENA
% measurements. Pairs with HAYSTACp2_process_v5
% - v4 created by DP, date uncertain
% - v5 created by DP on 9/8/19
% - v6 created by DP on 9/19/19
% - v7 created by DP on 9/23/19. Includes Y-factor processing
% - v8 created by DP on 9/23/19. Frequency-dependent correct for 1Q gains
% - v9 created by DP on 9/23/19. Relies entirely on function call for abg processing
% - v10 created by DP on 11/22/19. Takes into account the finitely detuned cavity for occasional
% calibration measurements as in ref. [4]
% - v11 created by DP on 12/2/19
% - v12 created by DP on 12/3/19. Begins work on the value-function approach to optimizing
% calibration parameters
% - v13 created by DP on 12/4/19. Switches value function to taking in over gain ratio
% - v14 created by DP on 12/21/19. Looks at h/c gain ratio as a second paramter of interest
% - v15 created by DP on 1/22/20. For use comparing with Kelly. Gets the code more in line with hers
% - v15a created by DP on 1/30/20. For use on 9/19/19 data and earlier. Special changes flagged with
% "CHG15a"
% - v15b created by DP on 2/16/20. For use on 9/11/19 data and earlier. Special changes flagged with
% "CHG15b"
% 
%% REFERENCES
% [1]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > HAYSTACp2 Processing Running
% Jounral > Processing Inputs"
% [2]: MM's lab NB: "measurements > scan rate enhancment > cold cavity spectro," top of page
% [3]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Squeezing/Hot Rod Calibration >
% Measurmenet Instructions for Kelly Backes"
% [4]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Finitely Detuned Cavity & 
% Calibration Measurements"
% [5]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Relating the Gains"
%

%% ============================================================================================== %%
%% INITIALIZATION
clear;
close all;
delete(findall(0,'type','figure','tag','TMWWaitbar'));
dt_start = datetime; % record start time
Settings_and_Defaults_v4;

%% ============================================================================================== %%
%% GLOBAL
global chi_mm2_fcn
chi_mm2_fcn = @(kap_dif, kap_sum, det) (kap_dif^2 + 4*det.^2)./(kap_sum^2 + 4*det.^2);

%% ============================================================================================== %%
%% USER INPUT
ver = 'v15b'; % should match script name. Not for the data to be loaded
data_run = '20190903';
dv_main = '1'; % data version

% booleans:
pt_monitor_for_add_noise = 0; % as of ref. [4], 11/23/19, I think it might be best not to pt_monitor
pt_exists_yfac = 01; % CHG15a - actually, this was in error: the pt always exists
print_timing = 01;
save_figs = 01;


% HAYSTAC universal parameters (same for all raw spectra, known externally):
UNIV.B0_T = 8; % magentic field
UNIV.V_m3 = 1.545 * 1E-3; % volume of the cavity excluding the tuning rod
UNIV.SQ_cav_eta = db2pow(-2.0/2); % SQ -> cav power transmission - need to refine this
UNIV.cav_AMP_eta = db2pow(-2.0/2); % cav -> AMP power transmission - need to refine this
UNIV.T_rod_K = 0.225; % best estimate from Kelly B., given to me 8/30/19
UNIV.T_f_K = 0.061; % given to me by Kelly B., 8/30/19
UNIV.T_hot_K = 0.333; % given to me by Kelly B., 10/8/19
UNIV.VTS_SQ_eta = db2pow(-1.035); % given to me by Kelly B., 12/19/19
UNIV.det_norm_tone_MHz = 2; % CHG15a (from 10E-3) 

% frequency ranges of interest:
% range over which to integrate the abg calibration measurements: 
UNIV.f_center_abg_Hz = 200E3;
UNIV.f_width_abg_Hz = 2E3; % symmetric 

f_lo_spec_Hz = 45E3;
f_hi_spec_Hz = 0.9E6; % 1.9E6; % CHG15a - smaller range

y_fac_cav_det_MHz = 2; % CHG15a

y_fac_gain_tone_det_MHz = 0.2; 


% value-fucntion approach to obtaining amplifier gain scale. The gain estimate starts at that
% measured (average) before and after the run. The abg and y-factor measuremens are locked in ratio
% with it. Criteria for N_H, S_c (i.e. p), and G_s must be met for some percent of data points in 
% the window of interest. Subject to those criteria, the haloscope noise is maximized. See ref. [5]
VF.crit.fracmatch = 1; % minimum fraction of frequencies that must match our criteria
VF.crit.pmin = 1/6; % min allowed participation ratio
VF.crit.pmax = 1/2; % max allowed participation ratio
VF.crit.NHmin_1q = 5; % min allowed HEMT-referred added noise
VF.crit.Gs_min_1q_pow = 0.02; % this is a lower bound for 1/squeezer gain
VF.crit.Gs_max_1q_pow = 1;

VF.hot_mult_min = 0.9;
VF.hot_mult_max = 1.1; % nice to construct the max and min roughly symmetrically around 1

VF.GA_addnoise_scale_min_1q_dB = 20; 
VF.GA_addnoise_scale_max_1q_dB = 30.5; 
VF.GA_spec_fixed_1q_dB = pow2db(650); 

VF.n_hot_mults = 99; % this must be odd
VF.n_GAs = 101; % nice to have this be odd

%% ============================================================================================== %%
%% CHECKS
if mod(VF.n_hot_mults,2) == 0 % even
    error('number of hot gain multipliers must be odd');
end

%% ============================================================================================== %%
%% DERIVED QUANTITIES
f_lo_abg_Hz = UNIV.f_center_abg_Hz - UNIV.f_width_abg_Hz/2;
f_hi_abg_Hz = UNIV.f_center_abg_Hz + UNIV.f_width_abg_Hz/2;

% value-fucntion gain quantities: 
VF.GA_addnoise_scale_test_1q_dB = linspace(VF.GA_addnoise_scale_min_1q_dB, ...
    VF.GA_addnoise_scale_max_1q_dB, VF.n_GAs);
VF.GA_addnoise_scale_test_1q_pow = db2pow(VF.GA_addnoise_scale_test_1q_dB);
VF.GA_spec_fixed_1q_pow = db2pow(VF.GA_spec_fixed_1q_dB);
VF.hot_mults = linspace(VF.hot_mult_min, VF.hot_mult_max, VF.n_hot_mults);

%% ============================================================================================== %%
%% LOADING INFO
% directory suffixes:
meta_suf = 'par'; % meta-data (parameters)
tx1_suf = 'tx'; % first cavity transmission measurement suffix
tx2_suf = 'tx2'; % second cavity transmission measurement suffix
rf_suf = 'rfl'; % cavity transmission suffix
spec_suf = 'psa'; % raw spectrum suffix
G_AMP_suf = 'jpaamp'; % AMP gain suffix - not used at present
abg_suf = 'abg'; % alpha-beta-gamma calibration measurements
add_noise_suf = 'yfactor'; % added noise measuement

% form factor:
ff_dir = 'form_fac_data';
ff_file = '4to5GHz_corrected.xlsx';

% refit couplings: 
% Kelly re-fits the couplings with a slightly better algorithm to what gets used to generate the par
% files. She has given me files containing the preferred values
rc_dir = 'refit_couplings';

%% ============================================================================================== %%
%% OUTPUT DIRECTORY
output_folder = ['HAYSTACp2_raw_saved_', data_run, '_', ver];
% create directory if necessary:
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

fig_folder = ['HAYSTACp2_wkspc_construction_figs_', data_run, '_', ver];
% create directory if necessary:
if ~exist(fig_folder, 'dir')
    mkdir(fig_folder);
end

%% ============================================================================================== %%
%% LOAD PARAMETERS + AUXILARY DATA
% load summary (paramter - par) meta-data:
par_wkspc = [data_run '_', dv_main, '_0_', meta_suf, '.mat'];
par = load([data_run, f, par_wkspc]).par; % this struct contains all the meta-data
n_OPs = length(par); % extract the number of operating points (i.e. of raw spectra)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in the first of some of these files and do some simple extractions and checks:
temp_tx1 = load([data_run, f, data_run '_', dv_main, '_1_', tx1_suf, '.mat']); 
temp_spec = load([data_run, f, data_run '_', dv_main, '_1_', spec_suf, '.mat']); 
n_IF_fs = length(temp_spec.meanavgps.singlesided_freqaxis); % necessary for allocation
% assumes linear spacing:
raw_res_Hz = temp_spec.meanavgps.singlesided_freqaxis(3) - ...
    temp_spec.meanavgps.singlesided_freqaxis(2); 
UNIV.sub_per_raw = raw_res_Hz * temp_spec.nAcq * temp_spec.acqInfo.SegmentSize / ...
            temp_spec.acqInfo.SampleRate;
        
clear(varname(temp_tx1), varname(temp_spec));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% form factors: 
ff_dat = readmatrix([ff_dir, f, ff_file]);
ff_f_corr_GHz = ff_dat(:, 1); % thermal contraction-corrected form factor frequencies
ff_C010 = ff_dat(:, 2); % corresponding form factors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% refit couplings:
rc_1 = load([rc_dir, f, 'datebetaQinfo_batch1', '.mat']).datebetaQinfo;
rc_2 = load([rc_dir, f, 'datebetaQinfo_batch2', '.mat']).datebetaQinfo;
rc_3 = load([rc_dir, f, 'datebetaQinfo_batch3', '.mat']).datebetaQinfo;
rc.beta = [rc_1.beta, rc_2.beta, rc_3.beta];
rc.spec_date = [rc_1.spectrum_date, rc_2.spectrum_date, rc_3.spectrum_date];
rc.spec_idx = [rc_1.spectrum_it, rc_2.spectrum_it, rc_3.spectrum_it];
rc.Q = [rc_1.Q, rc_2.Q, rc_3.Q];

%% ============================================================================================== %%
%% CREATE/SAVE WORKSPACES
% the alpha-beta-gamma and y-factor calibration measurements do not run every time. We just use
% the data for the closest OP that has one, always prefering the higher-numbered measruement in
% the event of a tie (since that measurement will be closer in time, as the calibrations happen
% before the data is taken)
abg_exists = zeros(1, n_OPs+1); % the +1 is becasue the last cals are named as at a new OP
add_noise_exists = zeros(1, n_OPs+1);
for i = 1:n_OPs+1
    if exist([data_run, f, data_run '_', dv_main, '_', num2str(i), '_', abg_suf, '.mat'], 'file')
        abg_exists(i) = 1; 
    end 
    if exist([data_run, f, data_run '_', dv_main, '_', num2str(i), '_', add_noise_suf, '.mat'], ...
            'file')
        add_noise_exists(i) = 1; 
    end
end
% the abg and added noise measurements should all be at the same OPs. Check that that is so:
if abg_exists ~= add_noise_exists
    error('unexpected difference in calibration operating points');
else
    cals_exist = abg_exists;
    clear(varname(abg_exists), varname(add_noise_exists));
end
    
cal_closest = NaN(1, n_OPs);
dec = 0.1;

valid_inds = find(cals_exist) - 1;
n_cals = length(valid_inds); 
for i = 1:n_OPs
    % we must subtract 1 from the vector of existing calibration locations, as the calibrations 
    % actually occur at the end of the previous operating point to the one at which they are listed:
    [~, ind] = min(abs(i + dec - valid_inds)); % adding a small number lets us round up
    cal_closest(i) = valid_inds(ind);
end

% save a workspace and also (non-proprietary) data files for each raw spectrum (i.e. each tuning)
% set an order to go through the spectra
order = 1:n_OPs;
order = setdiff(order, valid_inds);
order = [valid_inds, order];

wb = waitbar(0, 'creating operating point data workspaces...');
counting_i = 1;
while counting_i <= length(order)
    i = order(counting_i);
    cal_valid = 1; % init 1
    
    is_full_cal_op = 0;
    if ismember(i, valid_inds)
        is_full_cal_op = 1;
    end
    
    % these workspaces apply at every operating point: 
    tx1 = load([data_run, f, data_run '_', dv_main, '_', num2str(i), '_', tx1_suf, '.mat']); 
    tx2 = load([data_run, f, data_run '_', dv_main, '_', num2str(i), '_', tx2_suf, '.mat']); 
    rf = load([data_run, f, data_run '_', dv_main, '_', num2str(i), '_', rf_suf, '.mat']); 
    spec = load([data_run, f, data_run '_', dv_main, '_', num2str(i), '_', spec_suf, '.mat']); 
    G_AMP = load([data_run, f, data_run '_', dv_main, '_', num2str(i), '_', G_AMP_suf, '.mat']); 
    
    % these ones are the same for all spectra
    OP_data.B0_T = UNIV.B0_T; 
    OP_data.V_m3 = UNIV.V_m3; 
    OP_data.SQ_cav_eta = UNIV.SQ_cav_eta;
    OP_data.cav_AMP_eta = UNIV.cav_AMP_eta;
    OP_data.T_rod_K = UNIV.T_rod_K; 
    OP_data.sub_per_raw = UNIV.sub_per_raw;
    OP_data.W_width_MHz = 1E-6 * UNIV.f_width_abg_Hz;
    OP_data.W_center_MHz = 1E-6 * UNIV.f_center_abg_Hz;
    
    OP_data.G_AMP_1q_pow = db2pow(VF.GA_spec_fixed_1q_dB);

    % selet relevant IF frequency range 
    all_aq_freqs_Hz = spec.meanavgps.singlesided_freqaxis;
    [~, ind_lo_IF] = min(abs(all_aq_freqs_Hz - f_lo_spec_Hz));
    [~, ind_hi_IF] = min(abs(all_aq_freqs_Hz - f_hi_spec_Hz));
    range_IF = ind_lo_IF:ind_hi_IF;
    OP_data.f_IF_MHz = 1E-6 * all_aq_freqs_Hz(range_IF).';
    full_f_MHz = 1E-6 * all_aq_freqs_Hz.'; % for pulling out probe tone heights
    
    % ENA transmission data
    OP_data.ENA_f_tx1_GHz = tx1.f_GHz_tx; 
    OP_data.ENA_f_tx2_GHz = tx2.f_GHz_tx2; 
    OP_data.ENA_f_rf_GHz = rf.f_GHz_rfl; 
    OP_data.ENA_f_G_AMP_GHz = G_AMP.f_GHz_jpaamp; 

    OP_data.ENA_Ssw_t1 = tx1.I_tx + 1i*tx1.Q_tx;
    OP_data.ENA_Ssw_t2 = tx2.I_tx2 + 1i*tx2.Q_tx2;
    OP_data.ENA_Sss = rf.I_rfl + 1i*rf.Q_rfl;
    
    OP_data.f_cav_t1_GHz = par(i).Cavity_freq_GHz_tx;
    OP_data.f_cav_t2_GHz = par(i).Cavity_freq_GHz_tx2;
    OP_data.f_cav_r_GHz = par(i).Cavity_freq_GHz_rfl;
    OP_data.f_cav_r2_GHz = par(i).Cavity_freq_GHz_rfl2;
    
    OP_data.direct_squeezing = par(i).squeezing;

    % measured AMP gain. Not necessarily the value we will use
    G_AMP_meas_1q_pow = (par(i).Amp_gain + par(i).Amp_gain2) / 2;
    % OP_data.G_AMP_1q_pow = (par(i).Amp_gain + par(i).Amp_gain2) / 2;
   
    % raw spectrum, either probe-tone adjusted or not. Better to go with not, given relevant
    % timescales and tone sizes
    [~, ind_tone_spec] = min(abs(full_f_MHz - UNIV.det_norm_tone_MHz));
    tone_heights = spec.meanavgps.pt_power_est_list;
    OP_data.tone_height = mean(tone_heights);
    OP_data.tone_height_std = (std(tone_heights));
    OP_data.raw_spec_IF = spec.meanavgps.singlesided_powerspecavg(range_IF); % unweighted
    OP_data.tone_height = spec.meanavgps.singlesided_powerspecavg(ind_tone_spec);
    full_raw_spec_data = spec.meanavgps.singlesided_powerspecavg; % unweighted
    % OP_data.raw_spec_IF = spec.meanavgps.singlesided_ptweightspec(range_IF); % weighted
    % full_raw_spec_data = spec.meanavgps.singlesided_ptweightspec; % weighted

    
     % calculate couplings. Remember, coupling_factor = kappa_ext / kappa_loss    
    % get beta from refit coupligns object:
    rc_date_inds = zeros(size(rc.spec_date)); % init zeros
    rc_spec_inds = zeros(size(rc.spec_date)); % init zeros, same size
    for r = 1:length(rc.spec_date)
        rc_date_inds(r) = strcmp(num2str(rc.spec_date(r)), data_run);
        rc_spec_inds(r) = (rc.spec_idx(r) == i);
    end
    betaQ_ind = find(rc_date_inds & rc_spec_inds); % will be the single index where both are true
    beta = rc.beta(betaQ_ind);
    Q_load = rc.Q(betaQ_ind);        
    
    OP_data.cav_kappa_ext_t1_kHz = (1E6*OP_data.f_cav_t1_GHz / Q_load) / (1 + 1/beta);
    OP_data.cav_kappa_loss_t1_kHz = OP_data.cav_kappa_ext_t1_kHz / beta;
    OP_data.cav_kappa_ext_t2_kHz = (1E6*OP_data.f_cav_t2_GHz / Q_load) / (1 + 1/beta);
    OP_data.cav_kappa_loss_t2_kHz = OP_data.cav_kappa_ext_t2_kHz / beta;
    kappa_s_t1_kHz = OP_data.cav_kappa_ext_t1_kHz + OP_data.cav_kappa_loss_t1_kHz;
    kappa_d_t1_kHz = OP_data.cav_kappa_ext_t1_kHz - OP_data.cav_kappa_loss_t1_kHz;
    kappa_s_t2_kHz = OP_data.cav_kappa_ext_t2_kHz + OP_data.cav_kappa_loss_t2_kHz;
    kappa_d_t2_kHz = OP_data.cav_kappa_ext_t2_kHz - OP_data.cav_kappa_loss_t2_kHz;
    kappa_s_avg_kHz = (kappa_s_t1_kHz + kappa_s_t2_kHz)/2;
    kappa_d_avg_kHz = (kappa_d_t1_kHz + kappa_d_t2_kHz)/2;
    
    OP_data.T_fridge_K = UNIV.T_f_K;
    
    % fit the bandwidth of the AMP from the gain Lorentzian
    % FLAG - see below: 
    % note that the gain curve as taken is on resonance, and contains the dip of the cavity. Ideally
    % there would be a second AMP-off measurement, but for now we can hopefully extract the
    % relevant information from the refelction curve.
    AMPoff_meas_pow = abs(rf.I_rfl + 1i*rf.Q_rfl).^2;
    % the range for the rfl data is smaller, so interplate there
    AMPon_meas_pow_big = abs(G_AMP.I_jpaamp + 1i*G_AMP.Q_jpaamp).^2;
    AMPon_meas_pow = interp1(G_AMP.f_GHz_jpaamp, AMPon_meas_pow_big, rf.f_GHz_rfl.');
    AMP_gain_pow = AMPon_meas_pow ./ AMPoff_meas_pow;
    auto_scale_gain_fit = 0; % do not auto-scale, since this is a normalized gain where 1 means 1

    gain_out = FUNC_smart_gain_v1(rf.f_GHz_rfl, AMP_gain_pow, auto_scale_gain_fit);
    OP_data.BW_AMP_MHz = gain_out.w_GHz*1E3;
    
      
    % the form factor applies to the whole resonance. We get it by interpolation from the
    % simulations performed by Berkeley
    f_cav_avg_GHz = (OP_data.f_cav_t1_GHz + OP_data.f_cav_t2_GHz) / 2;
    OP_data.C_010 = interp1(ff_f_corr_GHz, ff_C010, f_cav_avg_GHz);

    % analyze for this OP or import from another OP the clibration data:
    if is_full_cal_op % is an operating point where we perfom a full cal. These come first 
        % added noise calibration :
        % + 1 in file name because cals actually belong to the previous OP
        add_noise = load([data_run, f, data_run '_', dv_main, '_', num2str(i+1), '_', ...
            add_noise_suf, '.mat']);
        % order as [cold, hot]
        T_VTS_K = [UNIV.T_f_K, UNIV.T_hot_K];
        G_A_meas_1q_pow = [add_noise.cold_gain, add_noise.hot_gain];
        if pt_monitor_for_add_noise
            disp('using probe tone-adjusted spectra for added noise spectra');
            S_out_addnoise_norm_1q = [add_noise.noise_amp_on_cold_pt_monitor, ...
                add_noise.noise_amp_on_hot_pt_monitor];
        else
            disp('NOT using probe tone-monitored spectra for added noise spectra'); %#ok<*UNRCH>
            S_out_addnoise_norm_1q = [add_noise.noise_amp_on_cold, add_noise.noise_amp_on_hot];
        end
        
        if pt_exists_yfac % CHG15a
            % crucial: the pt_monitored values are NOT on their own gain normalized. The acquisition
            % code Danielle wrote for these multiplies back in the mean of all the gain tones. I
            % will take these gains out - not by taking the measured gains themselves out, but by
            % taking out the probe tone itself. This applies WHETHER OR NOT we are using the
            % pt_monitor values, in fact, since the non-pt_monitored ones will also need to be
            % gain-normalized for the upcoming function calls
            [~, ind_tone] = min(abs(add_noise.freqs - UNIV.det_norm_tone_MHz*1E6));
            if ind_tone_spec ~= ind_tone
                error('tone indicies mismatched');
                % else just use same index
            end
            gap = 2;
            w_side = 3;
            range_pre = (ind_tone - gap - w_side):(ind_tone - gap - 1);
            range_post = (ind_tone + gap + 1):(ind_tone + gap + w_side);
            pt_addnoise_out = S_out_addnoise_norm_1q(ind_tone, :);
            nf_addnoise_out = mean(S_out_addnoise_norm_1q([range_pre, range_post], :));
            pt_addnoise_pow = pt_addnoise_out - nf_addnoise_out;
            % to keep things simple, make it so we just change the hot gain:
            % do this by dividng out 1st element:
            pt_addnoise_pow = pt_addnoise_pow / pt_addnoise_pow(1); 
        else % no probe tone available in Y-factor measurement, use JPA-measured gains
            pt_addnoise_pow = [1, add_noise.hot_gain / add_noise.cold_gain];
        end
        
        
        % normalize by the best proxy we have for something proportional to the gain (it's only the
        % relative normaliztion between hot and cold that will matter): 
        S_out_addnoise_norm_1q = S_out_addnoise_norm_1q ./ pt_addnoise_pow; 
        VF.GA_addnoise_test_1q_pow = VF.GA_addnoise_scale_test_1q_pow.' * pt_addnoise_pow;
        
        % make into a 3D structures where the 3rd dimension is multiplying the hot gain by different
        % numbers: 
        S_out_addnoise_norm_1q = S_out_addnoise_norm_1q .* ones(1,1,VF.n_hot_mults);
        VF.GA_addnoise_test_1q_pow = VF.GA_addnoise_test_1q_pow .* ones(1,1,VF.n_hot_mults);
        for k = 1:VF.n_hot_mults
            S_out_addnoise_norm_1q(:,2,k) = S_out_addnoise_norm_1q(:,2,k) / VF.hot_mults(k);
            VF.GA_addnoise_test_1q_pow(:,2,k) = VF.GA_addnoise_test_1q_pow(:,2,k) * VF.hot_mults(k); 
        end
        
        % now we want to construct the gains we will use for the abg measurements: a is the same as
        % the cold gains, b and g the same as each other, and as the run gain
        VF.G_abg_test_1q_pow = NaN(VF.n_GAs,3); % order of second index is 1~a, 2~b, 3~g
        VF.G_abg_test_1q_pow(:,1) = VF.GA_addnoise_test_1q_pow(:,1); % a same as cold
        VF.G_abg_test_1q_pow(:,2) = VF.GA_spec_fixed_1q_pow; % b same as run
        VF.G_abg_test_1q_pow(:,3) = VF.G_abg_test_1q_pow(:,2); % g same as b
        
        % this part is just for curioisity. I don't think that we are getting as good information on
        % the probe tone by looking at the values we saved for it (which are its bin summed with
        % adjacent ones, with no attempt at subtraction. What I do above is intelligently subtract 
        % some bins off to the side)
%         n_subspec_grps = length(add_noise.pts_for_moncold);
%         n_subspecs_per_grp = length(add_noise.pts_for_moncold{1});
%         pts_moncold = NaN(n_subspec_grps, n_subspecs_per_grp);
%         pts_monhot = NaN(n_subspec_grps, n_subspecs_per_grp);
%         for s = 1:n_subspec_grps
%             pts_moncold(s, :) = add_noise.pts_for_moncold{s};
%             pts_monhot(s, :) = add_noise.pts_for_monhot{s};
%         end
%         mean_pt_moncold = mean(pts_moncold, 'all');
%         mean_pt_monhot = mean(pts_monhot, 'all');
%         std_pt_moncold = std(pts_moncold, [], 'all');
%         std_pt_monhot = std(pts_monhot, [], 'all');
        
        
        % this calculation gives the added noise over the full range, but it should not be fully 
        % believed for the reason described in ref. [4]. What this will be used for is simply
        % calculating S_c over the first HWHM
        % step 1: start by getting N_H in a region where we cannot see the cavity: 
        
        % CHG15a - pass the AMP gain ENA measurement into the added noise analyzer 
        AMP_gain_meas.AMPon_f_GHz = G_AMP.f_GHz_jpaamp;
        AMP_gain_meas.AMPon_meas_pow_big = AMPon_meas_pow_big;
        S_c_known_1q_qta = 0; % set to 0 initially - a special value that sets S_c to S_f
        [OP_data.f4NH_IF_Hz, VF.N_Hp_test_1q_qta, ~] = FUNC_analyze_add_noise_v13b(OP_data, i, ...
            T_VTS_K, VF, S_c_known_1q_qta, y_fac_gain_tone_det_MHz, ...
            UNIV.det_norm_tone_MHz, y_fac_cav_det_MHz, UNIV.T_f_K, UNIV.VTS_SQ_eta, ...
            UNIV.SQ_cav_eta, UNIV.cav_AMP_eta, add_noise.freqs, S_out_addnoise_norm_1q, ...
            save_figs, fig_folder, AMP_gain_meas);
        % ^CHG15a - final argument added here
        % ^^CHG15b - 6 ENA arguments removed here
        
        % abg calibrations for hot road and squeezing:
        % + 1 in file name because cals actually belong to the previous OP
        abg = load([data_run, f, data_run '_', dv_main, '_', num2str(i+1), '_', abg_suf, '.mat']); 
        abg_gain_tone_det_MHz = y_fac_gain_tone_det_MHz; % same detuning for probe tone
   
        % CHG15a - add the run gain measurement to the abg cal since the ENA cals there don't have
        % AMP on: 
        abg.AMPon_f_GHz = G_AMP.f_GHz_jpaamp;
        abg.AMPon_meas_pow_big = AMPon_meas_pow_big;
        
        % step 2: calculate S_c over the entire range:
        % note that if plotted, G_s should really make no sense here, because we are not accounting
        % for the hot rod:
        all_plots_off = 01; % do not plot on this preliminary step 
        abg_test_out = FUNC_proc_abg_v5b(UNIV, OP_data, abg, i, abg_gain_tone_det_MHz, VF, ...
            save_figs, fig_folder, all_plots_off); 
      
        % step 3: now use the S_c calculated just above in order to re-caculate the added noise of
        % the HEMT over the full frequency range.
        S_c_known_1q_qta = abg_test_out.S_c_test_1q_qta; % set to test values just determined
        [OP_data.f4NH_IF_Hz, VF.N_Hp_test_1q_qta, ~] = FUNC_analyze_add_noise_v13b(OP_data, i, ...
            T_VTS_K, VF, S_c_known_1q_qta, y_fac_gain_tone_det_MHz, ...
            UNIV.det_norm_tone_MHz, y_fac_cav_det_MHz, UNIV.T_f_K, UNIV.VTS_SQ_eta, ...
            UNIV.SQ_cav_eta, UNIV.cav_AMP_eta, add_noise.freqs, S_out_addnoise_norm_1q, ...
            save_figs, fig_folder, AMP_gain_meas); 
        % ^CHG15a - final argument added here as well
        % ^^CHG15b - 6 ENA arguments removed here as well
        
        % step 4: with OP_data (incl. N_H) now updated, calcuate G_s (S_c should come out the same)
        all_plots_off = 01;   
        abg_test_out = FUNC_proc_abg_v5b(UNIV, OP_data, abg, i, abg_gain_tone_det_MHz, VF, ...
            save_figs, fig_folder, all_plots_off); 
        
        % include criteria for Sc
        S_f_1q_qta = FUNC_singlequad_S_of_T_f_v1(UNIV.T_f_K, f_cav_avg_GHz);
        S_r_1q_qta = FUNC_singlequad_S_of_T_f_v1(UNIV.T_rod_K, f_cav_avg_GHz);
        VF.crit.Scmin_1q_qta = VF.crit.pmin*(S_r_1q_qta - S_f_1q_qta) + S_f_1q_qta;
        VF.crit.Scmax_1q_qta = VF.crit.pmax*(S_r_1q_qta - S_f_1q_qta) + S_f_1q_qta;

        select_params = FUNC_valfct_v6(UNIV, OP_data, i, abg_test_out, VF, gain_out, save_figs, ...
            fig_folder);
        cal_valid = select_params.valid_point_found;
        
        if ~cal_valid
            % if cal invalid, remove it from list of valid cals and recalculate list of closest cals
            valid_inds(valid_inds == i) = []; % remove invalid cal
            n_cals = length(valid_inds);
            for m = 1:n_OPs
                [~, ind] = min(abs(m + dec - valid_inds)); % adding a small number lets us round up
                cal_closest(m) = valid_inds(ind);
            end
            order = [order, i]; %#ok<AGROW> % add cal back into line for analysis 
        end
        
        % fill in selected values:
        OP_data.N_Hp_1q_qta = select_params.N_Hp_1q_qta;
        OP_data.abg_out = select_params.abg_out; 
        
    else % not an operating point where we perform a full calibration
        % load the recently saved data for the nearest operating point where a cal was performed:       
        nearest_cal_OP = load([output_folder, f, 'dat_HAYSTACp2_OP', num2str(cal_closest(i))], ...
            varname(OP_data));
        % extract the relevant data from that loaded opearting point:
        OP_data.f4NH_IF_Hz = nearest_cal_OP.OP_data.f4NH_IF_Hz;
        OP_data.N_Hp_1q_qta = nearest_cal_OP.OP_data.N_Hp_1q_qta;
        OP_data.abg_out = nearest_cal_OP.OP_data.abg_out; 
    end
    
    if cal_valid 
        save([output_folder, f, 'dat_HAYSTACp2_OP', num2str(i)], varname(OP_data));
    end
    counting_i = counting_i + 1;
    waitbar(counting_i/n_OPs, wb);
end
close(wb);
disp([tab, 'workspaces created']);

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