%% DESCRIPTION
% a host of settings and defaults. Call this code at the beginning of a MATLAB program to set all
% these in place
%
%% HISTORY
% - v1 created by Dan Palken on 23 March 2018
% - v2 created by DP on 29 Oct 2018
% - v3 created by DP on 22 Apr 2019
%
%% NOTES
% to see the factory defaults, enter: get(groot, 'factory')
% to see the active defaults, enter: get(groot, 'defaults')
% to modify one, use: set(groot, 'default[SettingName]', [setting]);


%% ============================================================================================== %%
%% DOCK
set(0, 'DefaultFigureWindowStyle', 'docked');

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
neoncarrot = [255, 153, 51] / 255;
greenyellow = [153, 255, 51] / 255;	
mediumspringgreen = [51, 255, 153] / 255;
dodgerblue = [51, 153, 255] / 255;
wildstrawberry = [255, 51, 153] / 255;
gorse = [255, 255, 51] / 255;
redorange = [255, 51, 51] / 255;
saddlebrown = [139, 69, 19] / 255;
goldenrod = [218, 165, 32] / 255;
mediumseagreen = [60, 179, 113] / 255;
palegreen = [152, 251, 152] / 255;
lightgreen = [144, 238, 144] / 255;
indigo = [75, 0, 130] / 255;
lavender = [230, 230, 250] / 255;
violet = [238, 130, 238] / 255;

%% ============================================================================================== %%
%% SIZE
main_font_sz = 20; %32;
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
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1.0);
set(groot, 'defaultAxesTitleFontWeight', 'normal');
set(groot, 'defaultGeoaxesTitleFontWeight', 'normal');
set(groot, 'defaultPolaraxesTitleFontWeight', 'normal');
set(groot, 'defaultAxesFontWeight', 'normal');
set(groot, 'defaultUitableFontWeight', 'normal');
set(groot, 'defaultUipanelFontWeight', 'normal');
set(groot, 'defaultUicontrolFontWeight', 'normal');
set(groot, 'defaultUibuttongroupFontWeight', 'normal');
set(groot, 'defaultTextboxshapeFontWeight', 'normal');
set(groot, 'defaultTextarrowshapeFontWeight', 'normal');
set(groot, 'defaultTextFontWeight', 'normal');
set(groot, 'defaultColorbarFontWeight', 'normal');
set(groot, 'defaultGeoaxesFontWeight', 'normal');
set(groot, 'defaultGraphplotEdgeFontWeight', 'normal');
set(groot, 'defaultGraphplotNodeFontWeight', 'normal');
set(groot, 'defaultLegendFontWeight', 'normal');
set(groot, 'defaultPolaraxesFontWeight', 'normal');

%% #################################################################################################
%% ######################################## END OF PROGRAM #########################################
%% #################################################################################################