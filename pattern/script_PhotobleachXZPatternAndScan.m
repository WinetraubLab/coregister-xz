% Run this script to use Thorlabs system to scan a 3D OCT Volume and
% photobleach a given XZ pattern.
%
% This script performs:  
% 1. High-res OCT volume scan of tissue
% 2. Low-res surface detection scan
% 3. Precision XZ pattern photobleaching using the identified surface alignment  
%
% Configure the parameters in the INPUTS sections bellow:

%% INPUTS [1/2] - Perform Main OCT Scan

% Define the 3D Volume
pixel_size_um = 1; % x-y Pixel size in microns
xOverall_mm = [-0.25 0.25]; % Define the overall volume you would like to scan [start, finish]. OBJECTIVE_DEPENDENT: For 10x use [-0.5 0.5], for 40x use [-0.25 0.25]
yOverall_mm = [-0.25 0.25]; % Define the overall volume you would like to scan [start, finish]. OBJECTIVE_DEPENDENT: For 10x use [-0.5 0.5], for 40x use [-0.25 0.25]

% Define probe
octProbePath = yOCTGetProbeIniPath('40x','OCTP900'); % Inputs to the function are OBJECTIVE_DEPENDENT: '10x' or '40x', and scanning system dependent 'OCTP900' or ''
octProbeFOV_mm = 0.5; % How much of the field of view to use from the probe. OBJECTIVE_DEPENDENT: For 10x use 1, for 40x use 0.5
oct2stageXYAngleDeg = 0; % Angle between x axis of the motor and the Galvo's x axis

% Define z stack and z-stitching
scanZJump_um = 5; % microns. OBJECTIVE_DEPENDENT: For 10x use 15, for 40x use 5
zToScan_mm = unique([-100 (-30:scanZJump_um:300), 0])*1e-3; %[mm]
focusSigma = 10; % When stitching along Z axis (multiple focus points), what is the size of each focus in z [pixels]. OBJECTIVE_DEPENDENT: for 10x use 20, for 40x use 20 or 1

% Other scanning parameters
tissueRefractiveIndex = 1.33; % Use either 1.33 or 1.4 depending on the results. Use 1.4 for brain.
dispersionQuadraticTerm=-1.465e+08;  % 10x, OCTP900. This input is OBJECTIVE_DEPENDENT

% Where to save scan files
output_folder = 'E:\4.16\test';

% Set to true if you would like to process existing scan rather than scan a new one.
skipScanning = false;

% For all B-Scans, this parameter defines the depth (Z, pixels) that the focus is located at.
% If set to NaN, yOCTFindFocusTilledScan will be executed to request user to select focus position.
focusPositionInImageZpix = NaN;


%% INPUTS [2/2] - Photobleach XZ Pattern

% Define scanning & photobleaching area
skipHardware = false; % If true, it does not physically photobleach.

% Hashtag Photobleach Configurations
exposure_mm_sec = 5; % mm/sec
nPasses = 4; % Keep as low as possible. If galvo gets stuck, increase number

% Define X/Y ranges for surface detection. These determine the area scanned to map tissue topography 
% for depth correction. Set ranges to fully cover the photobleaching pattern area 
x_Surface_Detection_Range_mm = [-1.25, 1.25]; % (e.g., [-1.25, 1.25] = 2.5mm width)
y_Surface_Detection_Range_mm = [-1.25, 1.25];

% Generate the Hashtag Photobleaching Pattern (for full area)
[x_start_mm, x_end_mm, ...
 y_start_mm, y_end_mm, ...
 z_mm] = generateXZPattern();


%% Execution [1/3] - Perform Main OCT Scan

% Initialize a log file
log_file = Logging(output_folder, "initialize");
fprintf("\n[1/3] Performing Main OCT Scan...\n");
Logging(log_file, "update", "[1/3] Performing Main OCT Scan");

% Check that sufficient ammount of gel is above the tissue for proper focus
if (min(zToScan_mm)) > -100e-3
    warning('Because we use gel above tissue to find focus position. It is important to have at least one of the z-stacks in the gel. Consider having the minimum zToScan_mm to be -100e-3[mm]')
end

% Perform the OCT Scan
volumeOutputFolder = [output_folder '/OCTVolume/'];
fprintf('%s Please adjust the OCT focus such that it is precisely at the intersection of the tissue and the coverslip.\n', datestr(datetime));

fprintf('%s Scanning Volume\n',datestr(datetime));
scanParameters = yOCTScanTile (...
    volumeOutputFolder, ...
    xOverall_mm, ...
    yOverall_mm, ...
    'octProbePath', octProbePath, ...
    'tissueRefractiveIndex', tissueRefractiveIndex, ...
    'octProbeFOV_mm', octProbeFOV_mm, ...
    'pixelSize_um', pixel_size_um, ...
    'xOffset',   0, ...
    'yOffset',   0, ... 
    'zDepths',   zToScan_mm, ... [mm]
    'oct2stageXYAngleDeg', oct2stageXYAngleDeg, ...
    'skipHardware',skipScanning, ...
    'v',true  ...
    );

% Find focus in the scan if necessary
if isnan(focusPositionInImageZpix)
    fprintf('%s Find focus position volume\n',datestr(datetime));
    focusPositionInImageZpix = yOCTFindFocusTilledScan(volumeOutputFolder,...
        'reconstructConfig',{'dispersionQuadraticTerm',dispersionQuadraticTerm},'verbose',true);
end

% Reconstruct the z-stack 3d volume
fprintf('%s Processing\n',datestr(datetime));
outputTiffFile = [output_folder '/Image.tiff'];
yOCTProcessTiledScan(...
    volumeOutputFolder, ... Input
    {outputTiffFile},... Save only Tiff file as folder will be generated after smoothing
    'focusPositionInImageZpix', focusPositionInImageZpix,... No Z scan filtering
    'focusSigma',focusSigma,...
    'dispersionQuadraticTerm',dispersionQuadraticTerm,... Use default
    'interpMethod','sinc5', ...
    'v',true);

% Log scan & reconstruct completion
Logging(log_file, "update", "[1/3] Main OCT Scan Completed");
Logging(log_file, "update", "OCT Scan completed. If needed, you may now remove the hard drive.");
fprintf("[1/3] OCT Scan completed successfully.\n");


%% Execution [2/3] - Perform Low-Resolution Surface Identification Scans

Logging(log_file, "update", "[2/3] Performing Low-Resolution Surface Identification Scans");

[surfacePosition_mm, surface_x_mm, surface_y_mm, ~] = yOCTScanAndFindTissueSurface( ...
    'xRange_mm', x_Surface_Detection_Range_mm, ...
    'yRange_mm', y_Surface_Detection_Range_mm, ...
    'octProbeFOV_mm', octProbeFOV_mm, ...
    'octProbePath', octProbePath, ...
    'saveSurfaceMap', true, ...
    'output_folder', output_folder, ...
    'pixel_size_um', 10, ...
    'focusPositionInImageZpix', focusPositionInImageZpix, ...
    'dispersionQuadraticTerm', dispersionQuadraticTerm, ...
    'v', false);
surfaceFilePath = fullfile(output_folder, 'surface_data.mat');

Logging(log_file, "update", "[2/3] Finished Surface Identification Scans & Surface Map File Saved");
fprintf("[2/3] Surface Identification scans completed successfully.\n")


%% Execution [3/3] Photobleach XZ Pattern

Logging(log_file, "update", "[3/3] Photobleaching XZ Pattern (Hashtag Lines)");

yOCTPhotobleachTile(...
    [x_start_mm; y_start_mm],...
    [x_end_mm; y_end_mm],...
    'octProbePath',octProbePath,...
    'exposure',exposure_mm_sec,...
    'nPasses',nPasses,...
    'skipHardware',skipHardware, ...
    'z',z_mm, ...
    'surfaceMapFile', surfaceFilePath, ...
    'maxLensFOV', 0.5, ... Artificially reduce FOV to make sure lines are not tilted during photobleach
    'enableZoneAccuracy',1e-3,'enableZoneAccuracy',1e-3,... These arguments are needed for photobleaching the dots, comment out if no dotes are photobleached
    'plotPattern',true, ...
    'v',true);
disp('Done Photobleaching XZ Pattern')

Logging(log_file, "update", "[3/3] Finished Photobleaching XZ Pattern (Hashtag Lines)");
Logging(log_file, "update", "[3/3] XZ Scan & Pattern Script Ended Successfully!");


% Local function to log the process:
function log_filename = Logging(output_folder, mode, varargin)
    % Function to handle logging
    % - Mode 1: "initialize" -> Create log file and save the script contents
    % - Mode 2: "update" -> Append new process information to log
    % - Keeps using the same log file for all updates
    % - Accepts optional process_name for "update" mode

    persistent log_filepath; % Store log file path persistently

    % Default mode: If mode is not provided, set to "initialize"
    if nargin < 2 || isempty(mode)
        mode = "initialize";
    end

    % Get the path of the currently running script
    script_filename = matlab.desktop.editor.getActiveFilename(); 

    if isempty(script_filename)
        error('Could not retrieve the currently running script path.');
    end

    % Extract the directory where the script is running
    script_folder = fileparts(script_filename);

    % Default: If no output_folder is provided, use the script’s directory
    if nargin < 1 || isempty(output_folder)
        output_folder = script_folder;
    end

    % If initializing, create the log file and store the filename
    if strcmp(mode, "initialize")
        % Generate log file name
        log_index = 1;
        log_filepath = fullfile(output_folder, sprintf('runlog_%d.txt', log_index));

        % Ensure unique log file name
        while exist(log_filepath, 'file')
            log_index = log_index + 1;
            log_filepath = fullfile(output_folder, sprintf('runlog_%d.txt', log_index));
        end

        % Read the script’s contents
        script_text = fileread(script_filename);

        % Create log file and write script contents
        fid = fopen(log_filepath, 'w');
        fprintf(fid, "===== LOG FILE FOR RUN =====\n\n");
        fprintf(fid, "SCRIPT CONTENTS (FROM: %s):\n\n", script_filename);
        fprintf(fid, "%s\n\n", script_text);
        fclose(fid);

        fprintf("Log file created: %s\n", log_filepath);
        log_filename = log_filepath; % Return log filename

    elseif strcmp(mode, "update")
        % Ensure log_filepath is stored
        if isempty(log_filepath)
            error('Logging update called before initialize. Please run Logging in "initialize" mode first.');
        end

        % Default process_name if not provided
        if nargin < 3 || isempty(varargin)
            process_name = sprintf("[Process Update] Step Execution");
        else
            process_name = varargin{1};
        end

        % Append process details to the log
        fid = fopen(log_filepath, 'a');
        fprintf(fid, "---------------------------------------------------\n");
        fprintf(fid, "[%s] %s\n", datestr(datetime), process_name);
        fprintf(fid, "---------------------------------------------------\n\n");
        fclose(fid);

        fprintf("Log updated in: %s\n", log_filepath);
        log_filename = log_filepath; % Return log filename for consistency
    else
        error('Invalid mode. Use "initialize" or "update".');
    end
end
