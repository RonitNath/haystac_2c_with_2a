%% DESCRIPTION
% a host of settings and defaults. Call this code at the beginning of a MATLAB program to set all
% these in place
%
%% HISTORY
% - v1 created by Dan Palken on 23 March 2018
%
%% NOTES
% to see the factor defaults, enter: get(groot, 'factory')
% to modify one, use: set(groot, 'default[SettingName]', [setting]);

%% ============================================================================================== %%
%% COLOR
myred = [192, 80, 77] / 255;
myblue = [74, 126, 187] / 255;
mygreen = [119, 147, 60] / 255;
mypurple = [128, 100, 162] / 255;
mypink = [255, 47, 146] / 255;
myorange = [247, 150, 70] / 255;
mygray = [0.5, 0.5, 0.5];
navyblue = [0, 0, 0.8];
blue = [0.25, 0.25, 0.9];
red = [1, 0, 0];
black = [0, 0, 0];
white = [1, 1, 1];
yellow = [1, 1, 0];

%% ============================================================================================== %%
%% SIZE
main_font_sz = 22;
main_marker_sz = 15;

%% POSITION
default_left = 775;
default_bot = 25;
default_width = 650;
default_height = 475;
default_pos = [default_left, default_bot, default_width, default_height];

%% ============================================================================================== %%
%% SETTINGS
set(groot, 'defaultAxesFontName', 'latex');
set(groot, 'defaultTextFontName', 'latex');
set(groot, 'defaultRootFixedWidthFontName', 'latex');
set(groot, 'defaultTextarrowshapeFontName', 'latex');
set(groot, 'defaultPolaraxesFontName', 'latex');
set(groot, 'defaultTextboxshapeFontName', 'latex');
set(groot, 'defaultUibuttongroupFontName', 'latex');
set(groot, 'defaultUicontrolFontName', 'latex');
set(groot, 'defaultUipanelFontName', 'latex');
set(groot, 'defaultUitableFontName', 'latex');
set(groot, 'defaultColorbarFontName', 'latex');
set(groot, 'defaultLegendFontName', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');
set(groot, 'defaultPolaraxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextarrowshapeInterpreter', 'latex');
set(groot, 'defaultTextboxshapeInterpreter', 'latex');
set(groot, 'defaultAxesColorOrder', [black; black]);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesLineWidth', 4);
set(groot, 'defaultTextFontSize', main_font_sz);
set(groot, 'defaultAxesFontSize', main_font_sz);
set(groot, 'defaultColorbarFontSize', main_font_sz);
set(groot, 'defaultLegendFontSize', main_font_sz);
set(groot, 'defaultPolaraxesFontSize', main_font_sz);
set(groot, 'defaultTextboxshapeFontSize', main_font_sz);
set(groot, 'defaultAxesTickLength', 2.25 * [0.0100, 0.0250]);
set(groot, 'defaultAxesFontWeight', 'bold');
set(groot, 'defaultAxesBox', 'on');
set(groot, 'defaultErrorbarMarkerSize', main_marker_sz);
set(groot, 'defaultFunctionlineMarkerSize', main_marker_sz);
set(groot, 'defaultGraphplotMarkerSize', main_marker_sz);
set(groot, 'defaultImplicitfunctionlineMarkerSize', main_marker_sz);
set(groot, 'defaultLineMarkerSize', main_marker_sz);
set(groot, 'defaultAnimatedlineMarkerSize', main_marker_sz);
set(groot, 'defaultFunctionsurfaceMarkerSize', main_marker_sz);
set(groot, 'defaultImplicitfunctionsurfaceMarkerSize', main_marker_sz);
set(groot, 'defaultParameterizedfunctionlineMarkerSize', main_marker_sz);
set(groot, 'defaultParameterizedfunctionsurfaceMarkerSize', main_marker_sz);
set(groot, 'defaultPatchMarkerSize', main_marker_sz);
set(groot, 'defaultQuiverMarkerSize', main_marker_sz);
set(groot, 'defaultStairMarkerSize', main_marker_sz);
set(groot, 'defaultStemMarkerSize', main_marker_sz);
set(groot, 'defaultSurfaceMarkerSize', main_marker_sz);
set(groot, 'defaultSurfaceLineStyle', 'none');
set(groot, 'defaultFigureColorMap', jet);
set(groot, 'defaultFigurePosition', default_pos);

%% #################################################################################################
%% ######################################## END OF PROGRAM #########################################
%% #################################################################################################