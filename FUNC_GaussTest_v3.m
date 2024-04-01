function [delta] = FUNC_GaussTest_v3(counts, expected_counts, counts_upon_2_Poisson)

%% DESCRIPTION
% performs the test for Gaussianity as described in ref. [1]. In the proper limit, delta is the
% number of Gaussian standard deviations from the expected chi-squared-like figure of merit random
% variable. delta between +/-1 (or maybe +/-2) is pretty expected. Lower than -2 is either a rare
% chance fluctuation of the data towards precisely what the model would predict. Above 2 is either a
% rare chance fluctuation of the data away from expectation, or inidicative of the expected Gaussian
% not being the true distribution. 

% Some limitations: 
% It is not clear if this test can be used when the expected
% Gaussian comes from a fit (though for a high enough number of bins it might not make much
% difference). Also, there must be enough bins for our random variables to obey the CLT, or delta
% will not quite have the Gaussian interpretation we would like. Nonetheless, standard devations are
% standard deviations, and it should still be a reasonable figure of merit. 
%
%% HISTORY
% -v1 created by Dan Palken on 9/6/19
% -v2 created by DP on 9/8/19. Uses the prefeered statistic from ref. [1]
% -v3 created by DP on 9/8/19. Accounts for the fact that our histograms may have their counts/2 be
% Poisson distributed
%
%% REFERENCES
% [1]: DP's lab NB: "Data Acquisition + Analysis > HAYSTAC Phase 2 > Statistical Testing for
% Gaussianity"
%

%% ============================================================================================== %%
%% MAIN
if counts_upon_2_Poisson
    if sum(mod(counts, 2) ~= 0) > 0
        error('Input suggests counts are duplicated, but not all are even numbers'); 
    end
    expected_counts = expected_counts/2;
    counts = counts/2;
end

lambda = expected_counts;
N = length(lambda);
N_thresh = 30;
p_max_thresh = 0.05;


% define poisson parameters n, p from the expected distribution
n = sum(lambda);
p = lambda/n;
if N < N_thresh
    warning('number of bins is small: CLT may be vioalted');
end 
if max(p) > p_max_thresh
    warning(['large Poisson p = ', num2str(max(p), 2), ...
        ': may be outside of Poisson limit and into binomial']);
end

y = counts - lambda;
y2 = y.^2;
tot_y2 = sum(y2);

% expected statistics of the chi2 terms:
mu_y = lambda; 
var_y = 2 * lambda.^2 + lambda; % this +lambda is key for low-lambda terms - see ref. [1]
tot_mu_y = sum(mu_y);
tot_var_y = sum(var_y); % assumes independance - I believe reasonable
tot_std_y = sqrt(tot_var_y);

% final result: 
delta = (tot_y2 - tot_mu_y) / tot_std_y;
% - delta should be a mean-0 random variable
% - it may even be unbiased, as we are using an expected value for the standard deviation
% - delta = +(-) 1 corresponds to the data being 1 sigma noisy (1 sigma noiseless), etc...
% - under ideal enough circumstances (summing sufficeintly many terms), we might expect delta to be 
% standard normally distributed
% - if it's not standard nomral, we can still make weaker statements about what delta means - since
% its deviations from the mean in units of standard deviation - using Chebyshev's theorem, which
% says that at least 1-1/k^2 of the results will lie within k standard deviations (k > 1).

end

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################