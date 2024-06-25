function [...
    x_start_mm, x_end_mm, ...
    y_start_mm, y_end_mm, ...
    z_mm] = generateXZPattern()
% This function generates the instructions (lines) for XY lines
% The lines are from (x_start_mm(i), y_start_mm(i)) (to x_end_mm(i), y_end_mm(i)) 
% at depth z
% INPUT: verbose (default: false). If set to true, function will output
% line information as well as image 

%% Define the pattern 
base = 90e-3; %base seperation [mm]
lineLength = 6; %[mm]
vLineBias = 500e-3; % [mm]
hLineBias = -650e-3; % [mm]
vLinePositions = [-1 0 1]; % Unitless, vLine positions as multiplication of base
hLinePositions = [-3 0 1 3]; % Unitless, hLine positions as multiplication of base

z_mm = 30e-3; % [mm], how deep to draw pattern compared to the surface.

%% Build the pattern

% Define lines from pattern - v
x_start_mm = vLinePositions*base+vLineBias;
x_end_mm   = vLinePositions*base+vLineBias;
y_start_mm = -lineLength/2*ones(size(vLinePositions));
y_end_mm   = +lineLength/2*ones(size(vLinePositions));

% Define lines from pattern - h
x_start_mm = [x_start_mm -lineLength/2*ones(size(hLinePositions))];
x_end_mm   = [x_end_mm   +lineLength/2*ones(size(hLinePositions))];
y_start_mm = [y_start_mm  hLinePositions*base+hLineBias];
y_end_mm   = [y_end_mm    hLinePositions*base+hLineBias];

z_mm = z_mm * ones(size(x_start_mm));