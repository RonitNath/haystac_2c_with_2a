function [lag_vec, autocors] = FUNC_autocorr_fcn_v1(input, lag_extent)
%% DESCRIPTION:
% takes in some input vector and calculates the autocorrelation function (see the 'Estimation'
% section of ref. [1]) for all lags up to the specified extent with the lag equaling the index of
% the output vector
%
%% HISTORY:
% - v1 created by Dan Palken on 28 Jun 2018
%
%% REFERENCES
% [1]: https://en.wikipedia.org/wiki/Autocorrelation
%

%% ============================================================================================== %%
%% DERIVED QUANTITIES
n = length(input);
lag_vec = 1:lag_extent;
% estimate the mean and variance from the data itself:
input_var = var(input);
cent_input = input - mean(input); % subtract off mean

%% ============================================================================================== %%
%% ALLOCATE
autocors = zeros(1, lag_extent);

%% ============================================================================================== %%
%% BODY
for k = lag_vec % k is the lag
    prod = cent_input(1:n - k) .* cent_input(k + 1:n);
    autocors(k) = 1 / ((n - k - sum(isnan(prod))) * input_var) * nansum(prod);
end

end
% ##################################################################################################
% ######################################### END OF FUNCTION ########################################
% ##################################################################################################