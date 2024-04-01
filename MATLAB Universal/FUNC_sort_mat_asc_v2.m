function [out_mat] = FUNC_sort_mat_asc_v2(in_mat, ind_rows)
% DESCRIPTION:
% sorts a matrix in ascending order by the elements in one of its rows
%
% HISTORY:
% Created by: Dan Palken, 6 Sep 2017
%%
%% INPUTS:
%%  1. in_mat: the input matrix, generally unsorted
%%  2. ind_rows: a vector of (at least 1) row indecies that the matrix is
%%     to be sorted by. All sorts are ascending, and the first sort is by
%%     the outer-most index, with ties broken by the next-outermost, etc...
%%
%% OUTUTS:
%%  1. out-mat: the output matrix, now sorted - of the same dimensions as
%%     the input matrix
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DERIVED QUANTITIES
% number of times the data will have to be sorted
n_sorts = length(ind_rows);
n_cols = size(in_mat, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTER SORT
ind_row = ind_rows(1);

% gather indecies to reaggrange in order of
[~, sort_inds] = sort(in_mat(ind_row, :), 2);

out_mat = zeros(size(in_mat));
% perform the sort with the indecies gathered above
for i = 1:n_cols
    out_mat(:, i) = in_mat(:, sort_inds(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INNER SORTS
for i = 2:n_sorts
    % make a copy of the output matrix as it currently exists:
    out_mat_temp = out_mat;
    
    ind_row = ind_rows(i); % row to arrange by
    ind_prev_rows = ind_rows(1:(i - 1)); % all previous rows sorted by
    n_prev_rows = length(ind_prev_rows);
    
    % make a sub-matrix consisting of just those columns that are tied
    % after sorting by all previous rows
    sub_mat_start_ind = 1;
    sub_mat_stop_ind = 1;
    while sub_mat_stop_ind < n_cols
        
        tied_in_all_prev_rows = 1; % init true
        while tied_in_all_prev_rows
            
            % check if there is a tie in all previous rows
            for j = 1:n_prev_rows % loop over all previous rows sorted by
                ind_prev_row = ind_prev_rows(j); % specify a previous row
                if out_mat(ind_prev_row, sub_mat_start_ind) ~= ...
                        out_mat(ind_prev_row, sub_mat_stop_ind + 1)
                    tied_in_all_prev_rows = 0; % change to false
                end
            end
            
            % % increment stop index if all previous rows were tied
            if tied_in_all_prev_rows
                sub_mat_stop_ind = sub_mat_stop_ind + 1;  
            end
            
            % if the new stop index is the last possible index, simply set
            % the loop to stop by saying there was not a tie in all
            % previous rows
            if sub_mat_stop_ind == n_cols
                tied_in_all_prev_rows = 0; % change to false  
            end
            
            % then, if all previous rows were tied, loop back
        end
        
        % construct the sub-matrix
        sub_mat = out_mat(:, sub_mat_start_ind:sub_mat_stop_ind); 
        
        % sort that sub-matrix:
        % gather indecies to reaggrange in order of
        [~, sort_inds] = sort(sub_mat(ind_row, :), 2);
        
        sorted_sub_mat = zeros(size(sub_mat));
        % perform the sort with the indecies gathered above
        for k = 1:size(sub_mat, 2)
            sorted_sub_mat(:, k) = sub_mat(:, sort_inds(k));
        end
        
        % re-insert the sorted sub-matrix into its old spot in the output
        % matrix:
        out_mat(:, sub_mat_start_ind:sub_mat_stop_ind) = sorted_sub_mat; 

        
        % go to the next unchecked index before reentering the loop and
        % finding/sorting a new sub-matrix
        sub_mat_start_ind = sub_mat_stop_ind + 1;
        sub_mat_stop_ind = sub_mat_start_ind;
    end
end

end

% #########################################################################
% ############################ END OF FUNCTION ############################
% #########################################################################