%% DESCRIPTION
% Processes many data runs by calling wkspc_construction many times, using a few global variables
%
%% HISTORY
% - v1 created by Dan Palken on 1/25/20. Goes with wkspc_construction_v16
%

%% ============================================================================================== %%
%% INITIALIZATION
clear;
close all;
delete(findall(0,'type','figure','tag','TMWWaitbar'));
dt_start = datetime; % record start time

%% ============================================================================================== %%
%% MAIN
global data_run dv_main % declare globals

% boolean:
print_timing = 01;

% set globals and run workspace construction program
% try data_run = '20221103'; dv_main = '0'; wkspc_construction_v16; disp('test complete'); end
% data_run = '20221104'; dv_main = '0'; wkspc_construction_v16; disp('test complete');
% data_run = '20221106'; dv_main = '0'; wkspc_construction_v16; disp('test complete');
% data_run = '20221102'; dv_main = '0'; wkspc_construction_v16; disp('test complete');
% data_run = '20221026'; dv_main = '0'; wkspc_construction_v16; disp('test complete');
data_run = '20221111'; dv_main = '0'; wkspc_construction_v16; disp('test complete');
data_run = '20221114'; dv_main = '0'; wkspc_construction_v16; disp('test complete');
% data_run = '20190903'; dv_main = '1'; wkspc_construction_v16b; disp('20190903 complete');
% data_run = '20190907'; dv_main = '1'; wkspc_construction_v16b; disp('20190907 complete');
% data_run = '20190911'; dv_main = '0'; wkspc_construction_v16b; disp('20190911 complete');
% data_run = '20190915'; dv_main = '0'; wkspc_construction_v16a; disp('20190915 complete');
% data_run = '20190917'; dv_main = '2'; wkspc_construction_v16a; disp('20190917 complete');
% data_run = '20190919'; dv_main = '0'; wkspc_construction_v16a; disp('20190919 complete');
% data_run = '20190925'; dv_main = '0'; wkspc_construction_v16; disp('20190925 complete');
% data_run = '20190930'; dv_main = '0'; wkspc_construction_v16; disp('20190930 complete');
% data_run = '20191002'; dv_main = '4'; wkspc_construction_v16; disp('20191002 complete');
% data_run = '20191003'; dv_main = '0'; wkspc_construction_v16; disp('20191003 complete');
% data_run = '20191008'; dv_main = '2'; wkspc_construction_v16; disp('20191008 complete');
% data_run = '20191011'; dv_main = '4'; wkspc_construction_v16; disp('20191011 complete');
% data_run = '20191014'; dv_main = '0'; wkspc_construction_v16; disp('20191014 complete');
% data_run = '20191019'; dv_main = '2'; wkspc_construction_v16; disp('20191019 complete');
% data_run = '20191023'; dv_main = '0'; wkspc_construction_v16; disp('20191023 complete');
% data_run = '20191028'; dv_main = '0'; wkspc_construction_v16; disp('20191028 complete');
% data_run = '20191102'; dv_main = '0'; wkspc_construction_v16; disp('20191102 complete');
% data_run = '20191107'; dv_main = '0'; wkspc_construction_v16; disp('20191107 complete');
% data_run = '20191111'; dv_main = '0'; wkspc_construction_v16; disp('20191111 complete');
% data_run = '20191125'; dv_main = '6'; wkspc_construction_v16; disp('20191125 complete');
% data_run = '20191130'; dv_main = '0'; wkspc_construction_v16; disp('20191130 complete');
% data_run = '20191205'; dv_main = '0'; wkspc_construction_v16; disp('20191205 complete');
% data_run = '20191209'; dv_main = '0'; wkspc_construction_v16; disp('20191209 complete');
% data_run = '20191214'; dv_main = '0'; wkspc_construction_v16; disp('20191214 complete');

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