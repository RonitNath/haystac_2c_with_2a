function [S_singlequad_qta] = FUNC_singlequad_S_of_T_f_v1(T_K, f_GHz)
%% DESCRIPTION:
% takes in temperature (in K) and frequency (in GHz) and gives the single-quadrature spectral
% density of a matched load in units of quanta at those parameter values. Power is assumed to be
% divided equally between quadratures
%
%% HISTORY:
% - v1 created by Dan Palken on 4 May 2018
%

%% ============================================================================================== %%
%% CONVERT
% use the double-quadrature converter first, then divide by 2
S_doublqaud_qta = FUNC_S_of_T_f_v3(T_K, f_GHz);
S_singlequad_qta = S_doublqaud_qta / 2; % power gets divided between quadratures

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################