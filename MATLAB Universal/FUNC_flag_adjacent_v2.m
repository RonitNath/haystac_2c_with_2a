function [out_flags, ind_flagged_bins_out] = FUNC_flag_adjacent_v2(in_flags, n_adj)
%% DESCRIPTION:
% takes in a vector of 0's and 1's, where the 1's are assumed to be flags. Then flags the n_adj
% elements to the either side of every contiguous set of flags, and returns the whole set of flags
% as output. The upshot is that any set of m contiguous flags will come out as a set of (2 * n_adj +
% m) contiguous flags
%
%% HISTORY:
% - v1 created by Dan Palken on 9 Feb 2018
% - v2 created by DP on 22 Mar 2018: updated to newest coding conventions. Also corrected an error
% where row vector inputs would cause nonsense outputs  
%

%% ============================================================================================== %%
%% DERIVED QUANTITIES
ind_flagged_bins_in = find(in_flags); % locations of the flags
% crucially, this does not update as new flags are added
n_flags_orig = length(ind_flagged_bins_in);

%% ============================================================================================== %%
%% INITIALIZE
out_flags = in_flags; % start with only the input flags

%% ============================================================================================== %%
%% BODY
for i = 1:n_flags_orig
    f_ind = ind_flagged_bins_in(i); % index of the flag to be operated around
    
    left_surroundng_ind = f_ind - n_adj;
    right_surrounding_ind = f_ind + n_adj;
    
    % edge case: flag near left edge
    if f_ind <= n_adj
        left_surroundng_ind = 1;
    end
    
    % edge case: flag near right edge
    if f_ind > length(in_flags) - n_adj
        right_surrounding_ind = length(in_flags);
    end
    
   
    % plant flags around the original one:
    for j = left_surroundng_ind:right_surrounding_ind
        out_flags(j) = 1;
    end
end

ind_flagged_bins_out = find(out_flags);
end
% ##################################################################################################
% ######################################### END OF FUNCTION ########################################
% ##################################################################################################