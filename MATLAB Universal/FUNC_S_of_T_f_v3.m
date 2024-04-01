function [S_qta] = FUNC_S_of_T_f_v3(T_K, f_GHz)
%% DESCRIPTION:
% takes in temperature (in K) and frequency (in GHz) and gives the spectral density of a matched
% load in units of quanta at those parameter values
%
%% HISTORY:
% - v1 created by Dan Palken on 25 Aug 2017
% - v2 created by DP on 25 Sep 2017. T now in K
% - v3 created by DP on 19 Mar 2018. Updated to newest code conventions
%

%% ============================================================================================== %%
%% CONSTANTS
% constants of nature:
h_Js = 6.626070040E-34; % Planck constant
k_B_JperK = 1.38064852E-23; % Boltzmann constant

%% ============================================================================================== %%
%% DERIVED QUANTITIES
% put f in SI units: 
f_Hz = f_GHz * 1E9; 

%% ============================================================================================== %%
%% EQUATION
% note that the coth formula given here is a more compact version of the more intuitive 1 / 2 + ...
% formula
S_qta = (1 / 2) * coth((h_Js * f_Hz) ./ (2 * k_B_JperK * T_K));

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################