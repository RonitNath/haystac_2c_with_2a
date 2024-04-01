function [lag_vec, autocors] = FUNC_autocorr_fcn_v2(input, lag_extent)
%% DESCRIPTION:
% takes in some input vector and calculates the autocorrelation function (see the 'Estimation'
% section of ref. [1]) for all lags up to the specified extent with the lag equaling the index of
% the output vector
%
%% HISTORY:
% - v1 created by Dan Palken on 28 Jun 2018
% - v2 created by Dan Palken on 13 Sep 2019. Allows nans
%
%% REFERENCES
% [1]: https://en.wikipedia.org/wiki/Autocorrelation
%

%% ============================================================================================== %%
%% DERIVED QUANTITIES
n = length(input);
lag_vec = 1:lag_extent;
% estimate the mean and variance from the data itself:
input_var = nanvar(input);
cent_input = input - nanmean(input); % subtract off mean

%% ============================================================================================== %%
%% ALLOCATE
autocors = zeros(1, lag_extent);

%% ============================================================================================== %%
%% BODY
for k = lag_vec % k is the lag
    prod = cent_input(1:n-k) .* cent_input(1+k:n);
    n_nans_prod = sum(isnan(prod));
    autocors(k) = nansum(prod) / ((n-k-n_nans_prod) * input_var);
end

end
% ##################################################################################################
% ######################################### END OF FUNCTION ########################################
% ##################################################################################################