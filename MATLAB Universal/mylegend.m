function h = mylegend(varargin) 
%% DESCRIPTION
% wrapper used to set defaults
% adapted from ref. [1]
%
%% HISTORY
% - created by Dan Palken on 24 March 2018
%
%% REFERENCES
% [1]: https://www.mathworks.com/matlabcentral/answers/36640-default-text-size-in-legends

%% ============================================================================================== %%
%% USER INPUT
% default values:
lgd_font_sz = 18;
lgd_line_width = 3;

% boolean:
show_legends = 1; % option to disable showing legends

%% ============================================================================================== %%
%% BODY
if show_legends
    % declare object:
    h = legend(varargin{:});
    
    % set defaults:
    set(h, 'FontSize', lgd_font_sz);
    set(h, 'LineWidth', lgd_line_width);
else
    warning('mylegend disabled');
end

end

%% #################################################################################################
%% ######################################## END OF WRAPPER #########################################
%% #################################################################################################