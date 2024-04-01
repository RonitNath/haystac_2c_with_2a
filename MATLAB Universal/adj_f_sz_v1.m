function new_f_sz = adj_f_sz_v1(fig_width)

%% DESCRIPTION
% adjusts all font sizes to a new value with appropraite sizing for publications in mind (see ref.
% [1]. Input is the target width of the figure, and it adjusts all the sizes automatically. Returns
% the new font size
%
%% HISTORY
% - v1 created by Dan Palken on 18 October 2018
%
%% REFERENCES
% [1] DP's NB > Miscellaneous > Illustrator Tips 'n Tricks, 10/18/18 entry
%
%% NOTES
% to see the factor defaults, enter: get(groot, 'factory')
% to modify one, use: set(groot, 'default[SettingName]', [setting]);
%

%% ============================================================================================== %%
%% SIZE
default_size = 32; % 44;
default_width = 910;
new_f_sz = fig_width * (default_size / default_width); % calculated using ref. [1]

%% ============================================================================================== %%
%% ADJUST
set(groot, 'defaultTextFontSize', new_f_sz);
set(groot, 'defaultAxesFontSize', new_f_sz);
set(groot, 'defaultColorbarFontSize', new_f_sz);
set(groot, 'defaultLegendFontSize', new_f_sz);
set(groot, 'defaultPolaraxesFontSize', new_f_sz);
set(groot, 'defaultTextboxshapeFontSize', new_f_sz);

end
%% #################################################################################################
%% ######################################## END OF FUNCTION ########################################
%% #################################################################################################