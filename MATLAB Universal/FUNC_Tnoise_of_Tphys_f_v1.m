function [Tnoise_K] = FUNC_Tnoise_of_Tphys_f_v1(Tphys_K, f_GHz)
%% DESCRIPTION:
% takes in a physical temperature (in K) and frequency (in GHz) and gives the noise temperature,
% which is just the noise energy expressed in units of temperature
%
%% HISTORY:
% v1 created by Dan Palken on 11 Apr 2018 based somewhat on FUNC_S_of_T_f_v3
%
%% REFERENCES:
% [1]: DP's lab NB: "Data Acquisition + Analysis/Mock-Axion Experiment/Processing the
% Spectra/Processing the Spectra, part 2," 11 Apr 2018
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
% formula from ref. [3]
Tnoise_K = ((h_Js * f_Hz) / (2 * k_B_JperK)) * coth((h_Js * f_Hz) ./ (2 * k_B_JperK * Tphys_K));

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################