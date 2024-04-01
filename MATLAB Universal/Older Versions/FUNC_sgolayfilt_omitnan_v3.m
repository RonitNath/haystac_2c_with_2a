function output_spec = FUNC_sgolayfilt_omitnan_v3(input_spec, order, frame_len, freq_domain_fit)
%% DESCRIPTION:
% Applies a Savitsy-Golay filter function to allow a call where the input has NaN's, so that the
% NaN's are not included in the moving average of any of their peers. Output maintains all the
% positions of NaN's
%
%% HISTORY:
% - v1 created by Dan Palken on 9 Feb 2018. Based closely on Ling Zhong's mysgolayfit
% - v2 created by DP on 25 Mar 2018. Conventions updated and efficiency slightly improved
% - v3 created by DP on 25 Mar 2018. Option to apply filter in the time domain with fits at all NaN
% spots included. This should achieve nearly identical results to the frequency domain filtering,
% but orders of magnitude more quickly. Emperically, the two methods agree on each filtering to
% somewhere between a part in 500 and a part in 1,000 (that is, that is how far off the elements
% having the maximal discrepency are)
%

%% ============================================================================================== %%
% USER INPUT
% booleans:
% this code has the option to not put NaN's in the output where there were NaN's in the input. It
% will likely typically be preferable to keep those NaN's there, but to not include them, set to
% false:
output_NaNs = 1;
messages_on = 1;

% provide an update to the user at specified interval is messages are on
SG_update_percent = 4;

%% ============================================================================================== %%
%% DERIVED QUANTITIES
W = (frame_len - 1) / 2;
n_inputs = length(input_spec);

%% ============================================================================================== %%
%% ALLOCATE
output_spec = zeros(size(input_spec));

%% ============================================================================================== %%
% BODY
if messages_on
    disp('applying Savitzky-Golay filter...');
end

if freq_domain_fit % perform many, many polynomial fits
    SG_update = round(SG_update_percent * n_inputs / 100);
    
    for i = 1:n_inputs
        if i <= W % edge case: index near left edge
            x = (1:(i + W)).'; % indecies start at 1
            y = input_spec(1:(i + W));
        elseif i > n_inputs - W % edge case: index near right edge
            x = ((i - W):n_inputs).'; % indecies end at edge
            y = input_spec((i - W):n_inputs);
        else % main case: index in the middle somewhere
            x = ((i - W):(i + W)).';
            y = input_spec((i - W):(i + W)); % 2W + 1, centered on i
        end
        
        % perform the polynomial fit at the heart of the S-G filter, but remove all the terms where
        % there were NaN's. Since the x-values are kept, this does not effectively move some
        % elements adjacent to others that they should not be adjacent to
        p = polyfit(x(~isnan(y)), y(~isnan(y)), order);
        output_spec(i) = polyval(p, i); % i play the role of the x-value we evaluate at
        
        if mod(i, SG_update) == 0 && messages_on
            disp(['    SG filter: ', num2str(round(100 * i / n_inputs)), '% complete']);
        end
    end
    
else % no frequency domain fit: apply filter in the time domain 
    output_spec = input_spec; % initialzie as the input spectrum
    
    nan_locs = find(isnan(input_spec));
    n_nans = length(nan_locs);
    
    for j = 1:n_nans
        i = nan_locs(j); % get the index of the NaN which is being swapped out for a fitted value
        if i <= W % edge case: index near left edge
            x = (1:(i + W)).'; % indecies start at 1
            y = input_spec(1:(i + W));
        elseif i > n_inputs - W % edge case: index near right edge
            x = ((i - W):n_inputs).'; % indecies end at edge
            y = input_spec((i - W):n_inputs);
        else % main case: index in the middle somewhere
            x = ((i - W):(i + W)).';
            y = input_spec((i - W):(i + W)); % 2W + 1, centered on i
        end
        
        % perform the polynomial fit at the heart of the S-G filter, but only on the places where
        % there were NaN's in the original. Since this is a very time-costly thing to do, only doing
        % it a small percentage of the time is a big efficiency win
        p = polyfit(x(~isnan(y)), y(~isnan(y)), order);
        output_spec(i) = polyval(p, i); % i play the role of the x-value we evaluate at
    end
    
    % now take this output spectrum with all its NaN's repaced with fitted values and use the
    % standard MATLAB functionality
    output_spec = sgolayfilt(output_spec, order, frame_len);
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