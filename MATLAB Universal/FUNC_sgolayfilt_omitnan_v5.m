function [output_spec] = FUNC_sgolayfilt_omitnan_v5(input_spec, order, frame_len, freq_domain_fit)
%% DESCRIPTION:
% Applies a Savitsy-Golay filter function that allows a call where the input has NaN's. Output 
% maintains all the positions of NaN's
%
%% HISTORY:
% - v1 created by Dan Palken on 9 Feb 2018. Based closely on Ling Zhong's mysgolayfit
% - v2 created by DP on 25 Mar 2018. Conventions updated and efficiency slightly improved
% - v3 created by DP on 25 Mar 2018. Option to apply filter in the time domain with fits at all NaN
% spots included. This should achieve nearly identical results to the frequency domain filtering,
% but orders of magnitude more quickly. Emperically, the two methods agree on each filtering to
% somewhere between a part in 500 and a part in 1,000 (that is, that is how far off the elements
% having the maximal discrepency are)
% - v4 created by DP on 26 Mar 2018. Includes a thrid option which should never be used in practice,
% but is useful for quantifying the benefit on the intelligent time domain filtering implemented in
% v3. With this option, achieved by setting freq_domain_fit to 2, the NaN's are simply chucked out, 
% the standard sgolayfilt functionality is applied, and then the NaN's are reinserted. Emperically,
% this naive method gives errors of at most about 1 part in 200, so we can see that the more
% intelligent method implemented in v3 yields an improvement without sacrificing much on speed. The
% improvement seems to be about a factor of 4 in accuracy 
% - v5 created by DP on 29 Mar 2018. This version does a better job managing edge cases by always 
% calling my custom SG code (though there is a user option to call the native code), which takes a 
% tad longer than the native code, but produces more reasonable results near the ends of the
% spectrum. Also adds the ability to intelligently handle either row or column vector inputs, which 
% it could not do before for all modes of operation
%

%% ============================================================================================== %%
% USER INPUT
% booleans:
% this code has the option to not put NaN's in the output where there were NaN's in the input. It
% will likely typically be preferable to keep those NaN's there, but to not include them, set to
% false:
use_native_func = 1; % set to 0 to handle edge cases more intelligently
output_NaNs = 1;
% notifications:
messages_on = 0;
warnings_on = 0;

% provide an update to the user at specified interval is messages are on
SG_update_percent = 25;

%% ============================================================================================== %%
%% DERIVED QUANTITIES
W = (frame_len - 1) / 2;
n_inputs = length(input_spec);

%% ============================================================================================== %%
%% ALLOCATE
output_spec = zeros(size(input_spec));

%% ============================================================================================== %%
%% BODY
if messages_on
    disp('applying Savitzky-Golay filter...'); %#ok<*UNRCH>
end

if freq_domain_fit == 1 % perform many, many polynomial fits
    SG_update = round(SG_update_percent * n_inputs / 100);
    
    for i = 1:n_inputs
        if i <= W % edge case: index near left edge
            x = (1:(i + W)); % indecies start at 1
            y = input_spec(1:(i + W));
        elseif i > n_inputs - W % edge case: index near right edge
            x = ((i - W):n_inputs); % indecies end at edge
            y = input_spec((i - W):n_inputs);
        else % main case: index in the middle somewhere
            x = ((i - W):(i + W));
            y = input_spec((i - W):(i + W)); % 2W + 1, centered on i
        end
        
        if ~isequal(size(x), size(y)) % may need to transpose dimensions for the fit to work
            x = x.';
        end
        % perform the polynomial fit that the SG filter is equivalent to, but remove all the terms
        % where there were NaN's. Since the x-values are kept, this does not effectively move some
        % elements adjacent to others that they should not be adjacent to
        p = polyfit(x(~isnan(y)), y(~isnan(y)), order);
        output_spec(i) = polyval(p, i); % i plays the role of the x-value we evaluate at
        
        if mod(i, SG_update) == 0 && messages_on
            disp(['    SG filter: ', num2str(round(100 * i / n_inputs)), '% complete']);
        end
    end
    
elseif freq_domain_fit == 2 % do the naive thing and apply the standard functionality without NaN's
    if warnings_on
        % warn the user about this method: it is not intended for use in actual data processing:
        disp('    WARNING: ineffective Savitzky-Golay filter method selected');
    end
    
    nan_locs = find(isnan(input_spec));
    n_nans = length(nan_locs);

    output_spec = input_spec(~isnan(input_spec)); % remove the NaN's
    % apply the filter:
    if use_native_func
        if warnings_on
            % warn the user about poor edge processing:
            disp('    WARNING: suboptimal edge processing selected');
        end
        output_spec = sgolayfilt(output_spec, order, frame_len);
    else
        output_spec = FUNC_sgolayfilt_smartedge_v1(output_spec, order, frame_len);
    end
    % reinsert the NaN's:
    for j = 1:n_nans
        i = nan_locs(j); % index where the NaN was
        output_spec = [output_spec(1:i - 1), NaN, output_spec(i:end)];
    end
    
else % no frequency domain fit: apply filter in the time domain 
    output_spec = input_spec; % initialzie as the input spectrum
    
    nan_locs = find(isnan(input_spec));
    n_nans = length(nan_locs);
    
    for j = 1:n_nans
        i = nan_locs(j); % get the index of the NaN which is being swapped out for a fitted value
        if i <= W % edge case: index near left edge
            x = (1:(i + W)); % indecies start at 1
            y = input_spec(1:(i + W));
        elseif i > n_inputs - W % edge case: index near right edge
            x = ((i - W):n_inputs); % indecies end at edge
            y = input_spec((i - W):n_inputs);
        else % main case: index in the middle somewhere
            x = ((i - W):(i + W));
            y = input_spec((i - W):(i + W)); % 2W + 1, centered on i
        end
        
        if ~isequal(size(x), size(y)) % may need to transpose dimensions for the fit to work
            x = x.';
        end
        % perform the polynomial fit that the SG filter is equivalent to, but only on the places
        % where there were NaN's in the original. Since this is a very time-costly thing to do, only
        % doing it a small percentage of the time is a big efficiency win
        p = polyfit(x(~isnan(y)), y(~isnan(y)), order);
        output_spec(i) = polyval(p, i); % i plays the role of the x-value we evaluate at
    end
    
    % now take this output spectrum with all its NaN's repaced with fitted values and use the
    % my edge-case upgraded version of the native functionality (or the native functionality 
    % itself): 
    if use_native_func
        if warnings_on
            % warn the user about poor edge processing:
            disp('    WARNING: suboptimal edge processing selected');
        end
        output_spec = sgolayfilt(output_spec, order, frame_len);
    else 
        output_spec = FUNC_sgolayfilt_smartedge_v1(output_spec, order, frame_len);
    end
end

% optionally put NaNs in output where they originally were in input (generally this option will be
% kept on):
if output_NaNs
    output_spec(isnan(input_spec)) = NaN;
end

if messages_on
    disp('    Savitzky-Golay filter applied');
end

end

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################