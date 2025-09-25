%% photobleachSquareDermis.m - Square Pattern Photobleaching for Dermis Samples
%
% This script uses Thorlabs system to scan a 3D OCT Volume and photobleach
% a square XZ pattern with interior lines for dermis tissue samples.
%
% Workflow:
% 1. Low-res surface detection scan for tissue topography mapping
% 2. High-res OCT volume scan of tissue  
% 3. Precision XZ pattern photobleaching with surface alignment
% 4. Reconstruct 3D OCT volume from high-res Z-stack
%
% Prerequisites:
% - myOCT folder must be in MATLAB path
%   Example: addpath(genpath('F:\Jenkins\Scan OCTHist Dev\workspace\'))
% - Configure all parameters in the INPUTS section below

%% INPUTS

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
focusSigma = 10; % When stitching along Z axis (multiple focus points), what is the size of each focus in z [pixels]. OBJECTIVE_DEPENDENT: for 10x use 20, for 40x use 20 or 1

% Other scanning parameters
tissueRefractiveIndex = 1.33; % Use either 1.33 or 1.4 depending on the results. Use 1.4 for brain.

% Where to save scan files
outputFolder = '.\';

% Set to true if you would like to process existing scan rather than scan a new one.
skipHardware = true; % If true, skip real photobleaching and scanning

% How to apply Z correction from surface detection scan
surfaceCorrectionMode = 'none'; % Use 'per-tile' for best offset per tile/lens FOV, 'origin-tile' to use origin offset tile for all, or 'none' for no correction

%% Pre-processing and Pattern Configuration
volumeOutputFolder = [outputFolder '/OCTVolume/'];

% Photobleaching Configuration
nPasses = 4; % Keep as low as possible. If galvo gets stuck, increase number

% Pattern Parameters
L_sq_mm = 0.5;              % Square side length [mm]
pitch_mm = 0.02;            % Line spacing [mm] (20 µm)
z_list = (0:30:160)/1000;   % Z positions [mm] (0.03 to 0.18 mm)
centerX = 0.0;              % Pattern center X coordinate [mm]
centerY = 0.0;              % Pattern center Y coordinate [mm] 
includeBorders = true;      % Include border lines in square pattern

% Generate square pattern with interior lines
[x_start_mm, x_end_mm, y_start_mm, y_end_mm, z_mm] = generateXZSlabSingle( ...
    L_sq_mm, pitch_mm, z_list, centerX, centerY, 'v', includeBorders);

% Calculate pattern bounding box for surface detection
[xPhotobleachPatternBoundingBox, yPhotobleachPatternBoundingBox] = ...
    getPatternBoundingBox(x_start_mm, x_end_mm, y_start_mm, y_end_mm, octProbeFOV_mm);

%% Execution [0/4] - Focus and Dispersion Parameters

% Validate Z-scan range for proper focus detection
if (min(zToScan_mm)) > -100e-3
    warning(['Insufficient gel coverage detected. For proper focus detection, ', ...
             'at least one z-stack should be in the gel layer. ', ...
             'Consider setting minimum zToScan_mm to -100e-3 [mm]']);
end

fprintf('%s Please adjust OCT focus to tissue-coverslip intersection.\n', datestr(datetime));
fprintf('%s [0/4] Performing Focus and Dispersion Identification...\n', datestr(datetime));

% Focus and dispersion parameters (manually determined)
dispersionQuadraticTerm = -1.487424070517952e+08;
focusPositionInImageZpix = 447;

%% Execution [1/4] - Surface Detection Scan
fprintf('%s [1/4] Performing Low-Resolution Surface Detection...\n', datestr(datetime));

% Define main OCT scan bounding box [x y w h]
x0 = min(xOverall_mm);        % Left edge
y0 = min(yOverall_mm);        % Bottom edge  
w = abs(diff(xOverall_mm));   % Width
h = abs(diff(yOverall_mm));   % Height
mainScanBoundingBox_mm = [x0, y0, w, h];

% Skip surface scan when hardware is skipped or correction mode is 'none'
skipSurfaceScan = skipHardware || strcmpi(surfaceCorrectionMode, 'none');

% Perform tissue surface detection scan
[surfacePosition_mm, surfaceX_mm, surfaceY_mm] = yOCTTissueSurfaceAutofocus( ...
    'xRange_mm', xPhotobleachPatternBoundingBox, ...
    'yRange_mm', yPhotobleachPatternBoundingBox, ...
    'octProbeFOV_mm', octProbeFOV_mm, ...
    'octProbePath', octProbePath, ...
    'pixelSize_um', 10, ...
    'skipHardware', skipSurfaceScan, ...
    'focusPositionInImageZpix', focusPositionInImageZpix, ...
    'dispersionQuadraticTerm', dispersionQuadraticTerm, ...
    'roiToAssertFocus', mainScanBoundingBox_mm, ...
    'v', true);

% Store surface data in struct for easy passing
S.surfacePosition_mm = surfacePosition_mm;
S.surfaceX_mm = surfaceX_mm;
S.surfaceY_mm = surfaceY_mm;

%% Execution [2/4] - Main OCT Volume Scan
fprintf('%s [2/4] Performing Main OCT Volume Scan...\n', datestr(datetime));
fprintf('%s Scanning tissue volume...\n', datestr(datetime));
% (Optional) yOCTScanTile implementation here

%% Execution [3/4] - XZ Pattern Photobleaching
fprintf('%s [3/4] Starting XZ Pattern Photobleaching...\n', datestr(datetime));

% Pass 1: Square pattern with interior lines
fprintf('%s Photobleaching square pattern...\n', datestr(datetime));
exposure_sq = 3;  % Exposure time for square pattern
json_sq = yOCTPhotobleachTile( ...
    [x_start_mm; y_start_mm], ...
    [x_end_mm; y_end_mm], ...
    'octProbePath', octProbePath, ...
    'exposure', exposure_sq, ...
    'nPasses', nPasses, ...
    'skipHardware', skipHardware, ...
    'surfaceMap', S, ...
    'z', z_mm, ...
    'maxLensFOV', 0.4, ...
    'enableZoneAccuracy', 1e-3, ...
    'plotPattern', true, ...
    'v', true, ...
    'surfaceCorrectionMode', surfaceCorrectionMode);


% Pass 2: Thick cross pattern (confined to single FOV)
fprintf('%s Photobleaching thick cross pattern...\n', datestr(datetime));
maxLensFOV = 0.4;           % Maximum lens field of view [mm]
singleFOVMargin_mm = 0.01;  % Safety margin from FOV edges [mm]
L_cross_mm = min(L_sq_mm, maxLensFOV - 2*singleFOVMargin_mm); % Cross size <= 0.38 mm
thickOffset_mm = 0.01;      % Line thickness offset ~10 µm (3 lines per direction)

% Generate thick cross pattern
[xs2, xe2, ys2, ye2, z2] = generateXZCrossThickShort( ...
    L_cross_mm, z_list, centerX, centerY, thickOffset_mm);

exposure_cross = 6;  % Higher exposure for cross pattern
json2_cross = yOCTPhotobleachTile( ...
    [xs2; ys2], ...
    [xe2; ye2], ...
    'octProbePath', octProbePath, ...
    'exposure', exposure_cross, ...
    'nPasses', nPasses, ...
    'skipHardware', skipHardware, ...
    'surfaceMap', S, ...
    'z', z2, ...
    'maxLensFOV', maxLensFOV, ...
    'enableZoneAccuracy', 1e-3, ...
    'plotPattern', true, ...
    'v', true, ...
    'surfaceCorrectionMode', surfaceCorrectionMode);

%% Execution [4/4] - OCT Volume Reconstruction
fprintf('%s [4/4] Starting OCT Volume Reconstruction...\n', datestr(datetime));
outputTiffFile = [outputFolder '/Image.tiff'];
% (Optional) yOCTProcessTiledScan implementation here

%% Processing Complete
fprintf('%s Photobleaching and scanning completed successfully.\n', datestr(datetime));

%% Helper Functions

function [xBBox, yBBox] = getPatternBoundingBox(x_start_mm, x_end_mm, y_start_mm, y_end_mm, octProbeFOV_mm)
    % Calculate bounding box for photobleaching pattern
    % This determines the area scanned for tissue surface detection
    
    % Find pattern extents
    x_raw = [min([x_start_mm x_end_mm]), max([x_start_mm x_end_mm])];
    y_raw = [min([y_start_mm y_end_mm]), max([y_start_mm y_end_mm])];
    
    % Expand ranges to align with probe FOV
    expandRange = @(rng, FOV) [rng(1), rng(1) + ceil(diff(rng)/FOV) * FOV];
    
    xBBox = expandRange(x_raw, octProbeFOV_mm);
    yBBox = expandRange(y_raw, octProbeFOV_mm);
end

function [x_start_mm, x_end_mm, y_start_mm, y_end_mm, z_mm] = ...
    generateXZCrossThickShort(L_cross_mm, z_list_mm, centerX_mm, centerY_mm, thickOffset_mm)
    % Generate thick cross pattern confined to single FOV
    % Creates 3 parallel lines in each direction (vertical and horizontal)
    
    halfCross = L_cross_mm / 2;

    % Vertical lines (X fixed at 3 positions, Y spans cross length)
    xs_v = [centerX_mm - thickOffset_mm, centerX_mm, centerX_mm + thickOffset_mm];
    xv_start = xs_v;
    xv_end = xs_v;
    yv_start = (centerY_mm - halfCross) * ones(size(xs_v));
    yv_end = (centerY_mm + halfCross) * ones(size(xs_v));

    % Horizontal lines (Y fixed at 3 positions, X spans cross length)
    ys_h = [centerY_mm - thickOffset_mm, centerY_mm, centerY_mm + thickOffset_mm];
    xh_start = (centerX_mm - halfCross) * ones(size(ys_h));
    xh_end = (centerX_mm + halfCross) * ones(size(ys_h));
    yh_start = ys_h;
    yh_end = ys_h;

    % Combine vertical and horizontal lines
    x_start_mm = [xv_start, xh_start];
    x_end_mm = [xv_end, xh_end];
    y_start_mm = [yv_start, yh_start];
    y_end_mm = [yv_end, yh_end];

    % Assign Z positions cyclically across all lines
    nLines = numel(x_start_mm);
    if isscalar(z_list_mm)
        z_mm = z_list_mm * ones(1, nLines);
    else
        z_mm = zeros(1, nLines);
        for i = 1:nLines
            z_mm(i) = z_list_mm(1 + mod(i-1, numel(z_list_mm)));
        end
    end
end

