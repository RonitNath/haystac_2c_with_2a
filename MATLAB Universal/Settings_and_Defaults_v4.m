%% DESCRIPTION
% a host of settings and defaults. Call this code at the beginning of a MATLAB program to set all
% these in place
%
%% HISTORY
% - v1 created by Dan Palken on 23 March 2018
% - v2 created by DP on 29 Oct 2018
% - v3 created by DP on 22 Apr 2019
% - v4 created by DP on 20 Aug 2019. Uses objects
%
%% NOTES
% to see the factory defaults, enter: get(groot, 'factory')
% to see the active defaults, enter: get(groot, 'defaults')
% to modify one, use: set(groot, 'default[SettingName]', [setting]);

%% ============================================================================================== %%
%% SHORTCUT
f = filesep; % shortcut for platform-specific file separator
tab = '    '; % for dispaying messages

%% ============================================================================================== %%
%% DOCK OR NORMAL
% set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'DefaultFigureWindowStyle', 'normal');

%% ============================================================================================== %%
%% COLOR
SD.myred = [192, 80, 77] / 255;
SD.myblue = [74, 126, 187] / 255;
SD.mygreen = [119, 147, 60] / 255;
SD.mypurple = [128, 100, 162] / 255;
SD.myorange = [247, 150, 70] / 255;
SD.mypink = [255, 47, 146] / 255;
SD.mygray = [0.5, 0.5, 0.5];
SD.navyblue = [0, 0, 0.8];
SD.blue = [0.25, 0.25, 0.9];
SD.red = [1, 0, 0];
SD.black = [0, 0, 0];
SD.white = [1, 1, 1];
SD.yellow = [1, 1, 0];
SD.green = [0, 1, 0];
SD.orange = [1, 0.5, 0];
SD.purple = [0.5, 0, 0.5];
SD.neoncarrot = [255, 153, 51] / 255;
SD.greenyellow = [153, 255, 51] / 255;	
SD.mediumspringgreen = [51, 255, 153] / 255;
SD.dodgerblue = [51, 153, 255] / 255;
SD.wildstrawberry = [255, 51, 153] / 255;
SD.gorse = [255, 255, 51] / 255;
SD.redorange = [255, 51, 51] / 255;
SD.saddlebrown = [139, 69, 19] / 255;
SD.goldenrod = [218, 165, 32] / 255;
SD.mediumseagreen = [60, 179, 113] / 255;
SD.palegreen = [152, 251, 152] / 255;
SD.lightgreen = [144, 238, 144] / 255;
SD.indigo = [75, 0, 130] / 255;
SD.lavender = [230, 230, 250] / 255;
SD.violet = [238, 130, 238] / 255;

%% ============================================================================================== %%
%% SIZE
SD.main_font_sz = 17; % 20; %32;
SD.main_marker_sz = 15;

%% POSITION
SD.default_left = 775;
SD.default_bot = 25;
SD.default_width = 650;
SD.default_height = 475;
SD.default_pos = [SD.default_left, SD.default_bot, SD.default_width, SD.default_height];

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
set(groot, 'defaultAxesColorOrder', [SD.black; SD.black]);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesLineWidth', 4);
set(groot, 'defaultTextFontSize', SD.main_font_sz);
set(groot, 'defaultAxesFontSize', SD.main_font_sz);
set(groot, 'defaultColorbarFontSize', SD.main_font_sz);
set(groot, 'defaultLegendFontSize', SD.main_font_sz);
set(groot, 'defaultPolaraxesFontSize', SD.main_font_sz);
set(groot, 'defaultTextboxshapeFontSize', SD.main_font_sz);
set(groot, 'defaultAxesTickLength', 2.25 * [0.0100, 0.0250]);
set(groot, 'defaultAxesFontWeight', 'bold');
set(groot, 'defaultAxesBox', 'on');
set(groot, 'defaultErrorbarMarkerSize', SD.main_marker_sz);
set(groot, 'defaultFunctionlineMarkerSize', SD.main_marker_sz);
set(groot, 'defaultGraphplotMarkerSize', SD.main_marker_sz);
set(groot, 'defaultImplicitfunctionlineMarkerSize', SD.main_marker_sz);
set(groot, 'defaultLineMarkerSize', SD.main_marker_sz);
set(groot, 'defaultAnimatedlineMarkerSize', SD.main_marker_sz);
set(groot, 'defaultFunctionsurfaceMarkerSize', SD.main_marker_sz);
set(groot, 'defaultImplicitfunctionsurfaceMarkerSize', SD.main_marker_sz);
set(groot, 'defaultParameterizedfunctionlineMarkerSize', SD.main_marker_sz);
set(groot, 'defaultParameterizedfunctionsurfaceMarkerSize', SD.main_marker_sz);
set(groot, 'defaultPatchMarkerSize', SD.main_marker_sz);
set(groot, 'defaultQuiverMarkerSize', SD.main_marker_sz);
set(groot, 'defaultStairMarkerSize', SD.main_marker_sz);
set(groot, 'defaultStemMarkerSize', SD.main_marker_sz);
set(groot, 'defaultSurfaceMarkerSize', SD.main_marker_sz);
set(groot, 'defaultSurfaceLineStyle', 'none');
set(groot, 'defaultFigureColorMap', jet);
set(groot, 'defaultFigurePosition', SD.default_pos);
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