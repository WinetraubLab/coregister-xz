% Run this script to use Thorlabs system to scan a 3D OCT Volume and
% photobleach a given XZ pattern.
% This script performs:  
% 1. Low-res surface detection scan
% 2. High-res OCT volume scan of tissue
% 3. Precision XZ pattern photobleaching using the identified surface alignment
% 4. Reconstruct 3D OCT volume from high-res Z-stack
% Before running this script, make sure myOCT folder is in path for example
% by running: addpath(genpath('F:\Jenkins\Scan OCTHist Dev\workspace\'))
% Configure the parameters in the INPUTS sections bellow:

%% INPUTS
octSystem = 'Ganymede'; % Use either 'Ganymede' or 'Gan632' depending on your OCT system

% Define the 3D Volume
pixelSize_um = 1; % x-y Pixel size in microns
xOverall_mm = [-0.25 0.25]; % Define the overall volume you would like to scan [start, finish]. OBJECTIVE_DEPENDENT: For 10x use [-0.5 0.5], for 40x use [-0.25 0.25]
yOverall_mm = [-0.25 0.25]; % Define the overall volume you would like to scan [start, finish]. OBJECTIVE_DEPENDENT: For 10x use [-0.5 0.5], for 40x use [-0.25 0.25]

% Define probe
octProbePath = yOCTGetProbeIniPath('40x','OCTP900'); % Inputs to the function are OBJECTIVE_DEPENDENT: '10x' or '40x', and scanning system dependent 'OCTP900' or ''
octProbeFOV_mm = 0.5; % How much of the field of view to use from the probe. OBJECTIVE_DEPENDENT: For 10x use 1, for 40x use 0.5

% Define z stack and z-stitching
scanZJump_um = 5; % microns. OBJECTIVE_DEPENDENT: For 10x use 15, for 40x use 5
zToScan_mm = unique([-100 (-30:scanZJump_um:300), 0])*1e-3; %[mm]
cropZRange_mm = [min(zToScan_mm) max(zToScan_mm)]; % [mm] Crop output Z range: from above to below tissue surface. Set to [] to keep the full scan range without cropping.
focusSigma = 10; % When stitching along Z axis (multiple focus points), what is the size of each focus in z [pixels]. OBJECTIVE_DEPENDENT: for 10x use 20, for 40x use 20 or 1

% Other scanning parameters
tissueRefractiveIndex = 1.33; % Use either 1.33 or 1.4 depending on the results. Use 1.4 for brain.

% Hashtag Photobleach Configurations
exposure_mm_sec = 2; % mm/sec
nPasses = 4; % Keep as low as possible. If galvo gets stuck, increase number

% Where to save scan files
outputFolder = '.\';

% Set to true if you would like to process existing scan rather than scan a new one.
skipHardware = false; % If true, skip real photobleaching and scanning

% How to apply Z correction from surface detection scan
surfaceCorrectionMode = 'origin-tile'; % Use 'per-tile' for best offset per tile/lens FOV, 'origin-tile' to use origin offset tile for all, or 'none' for no correction

%% Load hardware
yOCTHardwareLibSetUp(octSystem, skipHardware, true);

%% Pre-processing
volumeOutputFolder = [outputFolder '/OCTVolume/'];

% Generate the Hashtag Photobleaching Pattern (for full area)
[x_start_mm, x_end_mm, y_start_mm, y_end_mm, z_mm] = generateXZPattern();

% Get the bounding box of the XZ pattern.
function [xPhotobleachPatternBoundingBox, yPhotobleachPatternBoundingBox] = ...
    getPatternBoundingBox(x_start_mm, x_end_mm, y_start_mm, y_end_mm, octProbeFOV_mm)

    % Define X/Y ranges for surface detection. These determine the area scanned to map tissue topography 
    % for depth correction. Set ranges to fully cover the photobleaching pattern area
    x_raw = [min([x_start_mm x_end_mm]), max([x_start_mm x_end_mm])]; % find X raw min and max
    y_raw = [min([y_start_mm y_end_mm]), max([y_start_mm y_end_mm])]; % find Y raw min and max
    expandRange = @(rng,FOV) [rng(1), rng(1) + ceil(diff(rng)/FOV) * FOV]; % helper that keeps the first edge, expands the second
    
    % Final areas for low-resolution surface scans
    xPhotobleachPatternBoundingBox = expandRange(x_raw, octProbeFOV_mm); % X Range e.g. [-1.25  1.25]
    yPhotobleachPatternBoundingBox = expandRange(y_raw, octProbeFOV_mm); % Y Range e.g. [-1.25  1.25]
end
[xPhotobleachPatternBoundingBox, yPhotobleachPatternBoundingBox] = ...
    getPatternBoundingBox(x_start_mm, x_end_mm, y_start_mm, y_end_mm, octProbeFOV_mm);

%% Execution [0/4] - Identify focus and dispersion parameters

% Check that sufficient ammount of gel is above the tissue for proper focus
if (min(zToScan_mm)) > -100e-3
    warning('Because we use gel above tissue to find focus position. It is important to have at least one of the z-stacks in the gel. Consider having the minimum zToScan_mm to be -100e-3[mm]')
end

fprintf('%s Please adjust the OCT focus such that it is precisely at the intersection of the tissue and the coverslip.\n', datestr(datetime));
fprintf('%s [0/4] Performing Focus and Dispersion Identification Scans...\n', datestr(datetime));

% Estimate dispersionQuadraticTerm and focusPositionInImageZpix using the
% glass slide.
if ~skipHardware
    [dispersionQuadraticTerm, focusPositionInImageZpix] = ...
        yOCTScanGlassSlideToFindFocusAndDispersionQuadraticTerm( ...
        'octProbePath', yOCTGetProbeIniPath('40x','OCTP900'), ...
        'tissueRefractiveIndex', tissueRefractiveIndex, ...
        'v', false);
    fprintf('%s dispersionQuadraticTerm = %.3e; focusPositionInImageZpix=%d;\n', ...
        datestr(datetime), dispersionQuadraticTerm, focusPositionInImageZpix);
end

% Uncomment below to set manually
% dispersionQuadraticTerm=-1.549e08;
% focusPositionInImageZpix = 200;

%% Execution [1/4] - Perform Low-Resolution Surface Identification Scans
fprintf('%s [1/4] Performing Low-Resolution Surface Identification Scans...\n', datestr(datetime));

% Define main OCT scan bounding box [x y w h]
x0 = min(xOverall_mm);        % left edge
y0 = min(yOverall_mm);        % bottom edge
w  = abs(diff(xOverall_mm));  % width
h  = abs(diff(yOverall_mm));  % height
mainScanBoundingBox_mm = [x0, y0, w, h];

% Run surface identification scan
if ~strcmpi(surfaceCorrectionMode,'none') % Skip surface scan when correction mode is 'none'
    [surfacePosition_mm, surfaceX_mm, surfaceY_mm] = yOCTTissueSurfaceAutofocus( ...
        'xRange_mm', xPhotobleachPatternBoundingBox, ...
        'yRange_mm', yPhotobleachPatternBoundingBox, ...
        'octProbeFOV_mm', octProbeFOV_mm, ...
        'octProbePath', octProbePath, ...
        'pixelSize_um', 10, ...
        'focusPositionInImageZpix', focusPositionInImageZpix, ...
        'dispersionQuadraticTerm', dispersionQuadraticTerm, ...
        'roiToAssertFocus', mainScanBoundingBox_mm, ...
        'v', true);
    S.surfacePosition_mm = surfacePosition_mm;
    S.surfaceX_mm        = surfaceX_mm;
    S.surfaceY_mm        = surfaceY_mm;
else
    S = [];
end

%% Execution [2/4] - Perform Main OCT Scan
% Perform the OCT Scan
fprintf('%s [2/4] Performing Main OCT Scan...\n', datestr(datetime));
fprintf('%s Scanning Volume\n',datestr(datetime));
scanParameters = yOCTScanTile (...
    volumeOutputFolder, ...
    xOverall_mm, ...
    yOverall_mm, ...
    'octProbePath', octProbePath, ...
    'tissueRefractiveIndex', tissueRefractiveIndex, ...
    'octProbeFOV_mm', octProbeFOV_mm, ...
    'pixelSize_um', pixelSize_um, ...
    'xOffset',   0, ...
    'yOffset',   0, ... 
    'zDepths',   zToScan_mm, ... [mm]
    'v',true  ...
    );

%% Execution [3/4] Photobleach XZ Pattern
fprintf('%s [3/4] Starting Photobleaching XZ Pattern\n', datestr(datetime));
json = yOCTPhotobleachTile(...
            [x_start_mm; y_start_mm],...
            [x_end_mm; y_end_mm],...
            'octProbePath',octProbePath,...
            'exposure',exposure_mm_sec,...
            'nPasses',nPasses,...
            'surfaceMap', S, ...
            'z',z_mm, ...
            'maxLensFOV', 0.4 , ... mm Artificially reduce FOV to make sure lines are not tilted during photobleach
            'enableZoneAccuracy',1e-3,'enableZoneAccuracy',1e-3,... These arguments are needed for photobleaching the dots, comment out if no dotes are photobleached
            'plotPattern',true, ...
            'v',true,...
            'surfaceCorrectionMode', surfaceCorrectionMode);

%% Cleanup for next run
yOCTHardwareLibTearDown(true);

%% Main OCT Volume Reconstruction [4/4]
% Reconstruct the z-stack 3d volume
fprintf('%s [4/4] Starting Reconstruction.\n', datestr(datetime));

outputTiffFile = [outputFolder '/Image.tiff'];
if ~skipHardware
    yOCTProcessTiledScan(...
        volumeOutputFolder, ... Input
        {outputTiffFile},... Save only Tiff file as folder will be generated after smoothing
        'focusPositionInImageZpix', focusPositionInImageZpix,... No Z scan filtering
        'focusSigma',focusSigma,...
        'cropZRange_mm',cropZRange_mm,...
        'outputFilePixelSize_um', pixelSize_um,...
        'dispersionQuadraticTerm',dispersionQuadraticTerm,... Use default
        'interpMethod','sinc5', ...
        'v',true);
end

%% Clean up
fprintf('%s All done.\n', datestr(datetime));
