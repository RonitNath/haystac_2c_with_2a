function output_spec = FUNC_sgolayfilt_omitnan_v2(input_spec, order, framelen)
%% DESCRIPTION:
% Applies a Savitsy-Golay filter function to allow a call where the input has NaN's, so that the
% NaN's are not included in the moving average of any of their peers. Output maintains all the
% positions of NaN's. 
%
% Note that the reason we cannot do the obvious thing and just acheive this using a weight vector as
% the standard MATLAB functionality seems to permit is that the weight vector it intakes is for the
% moving window, not for the input vector
%
%% HISTORY:
% - v1 created by Dan Palken on 9 Feb 2018. Based closely on Ling Zhong's mysgolayfit
% - v2 created by DP on 25 Mar 2018. Conventions updated and efficiency slightly improved
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
% DERIVED QUANTITIES
W = (framelen - 1) / 2;
n_inputs = length(input_spec);
SG_update = round(SG_update_percent * n_inputs / 100);

%% ============================================================================================== %%
% ALLOCATE
output_spec = zeros(size(input_spec));

%% ============================================================================================== %%
% BODY
if messages_on
    disp('applying Savitzky-Golay filter...');
end
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
    % there were NaN's. Since the x-values are kept, this does not effectively move some elements
    % adjacent to others that they should not be adjacent to
    p = polyfit(x(~isnan(y)), y(~isnan(y)), order);
    output_spec(i) = polyval(p, i);
    
    if mod(i, SG_update) == 0 && messages_on
        disp(['    SG filter: ', num2str(round(100 * i / n_inputs)), '% complete']);
    end
end

% optionally put NaNs in output where they originally were in input:
if output_NaNs
    output_spec(isnan(input_spec)) = NaN;
end

end

%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################