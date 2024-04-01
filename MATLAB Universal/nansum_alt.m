function [s, all_nan_locs, n_all_nan_locs] = nansum_alt(x)

%% DESCRIPTION
% performs a nansum, but if all terms are nan, instead of providing zeros, gives NaN. Tested for
% inputs up to 2 dimensions. No guarantees that it works higher than that, but it might
%
%% HISTORY
% - v1 created by Dan Palken on 13 Sep 2019. . For now, no bells and whistles
%

%% ============================================================================================== %%
%% MAIN
% this may or may not be necessary, but it's quick, it makes things easy to visualize, and it works
if size(x, 1) == 1 % make column vec, so we sum along first dimension
    x = x.';
end

s = nansum(x, 1); % sum along first dimension - if x is a matrix will get a row vector

nan_locs = isnan(x);
% find locations in the remaining structure where there were all NaN's that went into the sum:
all_nan_locs = sum(nan_locs, 1) == size(x, 1); 
s(all_nan_locs) = NaN; % set those locations to NaN's
n_all_nan_locs = sum(all_nan_locs);

end
%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################