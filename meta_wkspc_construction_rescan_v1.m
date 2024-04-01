%% DESCRIPTION
% Processes many data runs by calling wkspc_construction many times, using a few global variables
%
%% HISTORY
% - v1 created by Dan Palken on 3/24/20. Based riectly off of meta_wkspc_construction_v1. Goes with 
% wkspc_construction_rescan_v2
%

%% ============================================================================================== %%
%% INITIALIZATION
clear;
close all;
addpath('MATLAB Universal')
delete(findall(0,'type','figure','tag','TMWWaitbar'));
dt_start = datetime; % record start time

%% ============================================================================================== %%
%% MAIN
global data_run dv_main % declare globals

% boolean:
print_timing = 01;

% set globals and run workspace construction program
data_run = '20200221'; dv_main = '3'; wkspc_construction_rescan_v2; disp('20200221 complete');
data_run = '20200225'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200225 complete');
data_run = '20200226'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200226 complete');
data_run = '20200227'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200227 complete');
data_run = '20200228'; dv_main = '1'; wkspc_construction_rescan_v2; disp('20200228 complete');
data_run = '20200229'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200229 complete');
data_run = '20200301'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200301 complete');
data_run = '20200303'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200303 complete');
data_run = '20200304'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200304 complete');
data_run = '20200305'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200305 complete');
data_run = '20200306'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200306 complete');
data_run = '20200307'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200307 complete');
data_run = '20200308'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200308 complete');
data_run = '20200309'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200309 complete');
data_run = '20200310'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200310 complete');
data_run = '20200311'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200311 complete');
data_run = '20200312'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200312 complete');
data_run = '20200313'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200313 complete');
data_run = '20200314'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200314 complete');
data_run = '20200315'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200315 complete');
data_run = '20200316'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200316 complete');
data_run = '20200317'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200317 complete');
data_run = '20200318'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200318 complete');
data_run = '20200319'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200319 complete');
data_run = '20200320'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200320 complete');
data_run = '20200321'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200321 complete');
data_run = '20200322'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200322 complete');
data_run = '20200323'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200323 complete');
data_run = '20200325'; dv_main = '1'; wkspc_construction_rescan_v2; disp('20200325 complete');
data_run = '20200326'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200326 complete');
data_run = '20200328'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200328 complete');
data_run = '20200329'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200329 complete');
data_run = '20200330'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200330 complete');
data_run = '20200331'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200331 complete');
data_run = '20200401'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200401 complete');
data_run = '20200402'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200402 complete');
data_run = '20200403'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200403 complete');
data_run = '20200404'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200404 complete');
data_run = '20200405'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200405 complete');
data_run = '20200406'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200406 complete');
data_run = '20200407'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200407 complete');
data_run = '20200408'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200408 complete');
data_run = '20200410'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200410 complete');
data_run = '20200411'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200411 complete');
% these last two are not axion-sensitive. They are simply to see if we can average down to KSVZ
% levels under normal operating conditions if we wanted to
%data_run = '20200414'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200414 complete');
%data_run = '20200417'; dv_main = '0'; wkspc_construction_rescan_v2; disp('20200417 complete');


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