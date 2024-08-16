% Run this demo to use Thorlabs system to photobleach a pattern of a square
% with an L shape on its side

% Before running this script, make sure myOCT folder is in path for example
% by running: addpath(genpath('F:\Jenkins\Scan OCTHist Dev\workspace\'))

%% Inputs

% When set to true the stage will not move and we will not
% photobleach. Use "true" when you would like to see the output without
% physcaily running the test.
skipHardware = false;

% OCT probe
octProbePath = yOCTGetProbeIniPath('40x','OCTP900'); % Select lens magnification

% Pattern to photobleach. System will photobleach n lines from 
% (x_start(i), y_start(i)) to (x_end(i), y_end(i)) at height z
% Load the pattern
[x_start_mm, x_end_mm, y_start_mm, y_end_mm, z_mm] = ...
    generateXZPattern();

% Photobleach configurations
exposure_mm_sec = 5; % mm/sec
nPasses = 4; % Keep as low as possible. If galvo gets stuck, increase number

%% Photobleach
yOCTPhotobleachTile(...
    [x_start_mm; y_start_mm],...
    [x_end_mm; y_end_mm],...
    'octProbePath',octProbePath,...
    'exposure',exposure_mm_sec,...
    'nPasses',nPasses,...
    'skipHardware',skipHardware, ...
    'z',z_mm(1), ...
    'maxLensFOV', 0.15, ... Artificially reduce FOV to make sure lines are not tilted during photobleach
    'enableZoneAccuracy',1e-3,'enableZoneAccuracy',1e-3,... These arguments are needed for photobleaching the dots, comment out if no dotes are photobleached
    'plotPattern',true, ...
    'v',true); 
disp('Done Patterning')

%% OCT Volume Scan
error('OCT Volume Scan Not Yet Defined')
