function output_spec = FUNC_sgolayfilt_smartedge_v1(input_spec, order, frame_len)
%% DESCRIPTION:
% Applies a Savitsy-Golay filter function that uses MATLAB's native functionality for all points not
% within half a frame-length-minus-one of the edge, and then manually does the fit for the remaining
% edge cases, as what MATLAB naturally does is frankly not all that smart (see ref. [1]). This
% functionality assumes there are no NaN's in the input spectrum. Call this (and not the native
% MATLAB sgolayfilt with e.g. the FUNC_sgolayfilt_omitnan_v4 code in order to handle NaN's
% appropriately
% 
%% HISTORY:
% - v1 created by Dan Palken on 29 Mar 2018. Uses some logic similar to that in e.g. 
% FUNC_sgolayfilt_omitnan_v4
%
%% REFERENCES:
% [1]: DP's lab NB, 'Data Acquisition + Analysis/Mock-Axion Experiment/Coding Tips & Trick/MATLAB's
% Savitsky-Golay Functionality,' 3/29/18 entry 
%

%% ============================================================================================== %%
%% DERIVED QUANTITIES
W = (frame_len - 1) / 2;
n_inputs = length(input_spec);

%% ============================================================================================== %%
%% BODY
% apply the native functionality. All of the values that are not within W of the edges will be 
% correct:
output_spec = sgolayfilt(input_spec, order, frame_len);

% go back and write over the edge cases:
for i = 1:n_inputs
    perform_fit = 1; % init true
    if i <= W % edge case: index near left edge
        x = (1:(i + W)); % indecies start at 1
        y = input_spec(1:(i + W));
    elseif i > n_inputs - W % edge case: index near right edge
        x = ((i - W):n_inputs); % indecies end at edge
        y = input_spec((i - W):n_inputs);
    else  % for main cases, do not fit
        perform_fit = 0;
    end
    
    if perform_fit % perform the polynomial fit that the S-G filter is equivalent to
        if ~isequal(size(x), size(y)) % may need to transpose dimensions for the fit to work
            x = x.';
        end
        p = polyfit(x, y, order);
        output_spec(i) = polyval(p, i);
    end
end

end

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################