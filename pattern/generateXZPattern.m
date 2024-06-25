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
vLineBias = 500e-3; % [mm]
hLineBias = -650e-3; % [mm]
base = 90e-3; %base seperation [mm]
vLinePositions = [-1 0 1]; % Unitless, vLine positions as multiplication of base
hLinePositions = [-3 0 1 3]; % Unitless, hLine positions as multiplication of base
lineLength = 6; %[mm]

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

%% Add dots
number_of_dots = 10;
rng(13); % Seed
rand_numbers_m1to1 = (2*rand(2, number_of_dots)-1); % Random numbers between -1 and 1
rand_numbers_m1to1 = round(rand_numbers_m1to1*4)/4; % Only specific cells are allowed
dot_xy_mm = rand_numbers_m1to1*0.1;

x_start_mm = [x_start_mm (dot_xy_mm(1,:)-1e-3)];
x_end_mm = [x_end_mm (dot_xy_mm(1,:)+1e-3)];
y_start_mm = [y_start_mm dot_xy_mm(2,:)];
y_end_mm = [y_end_mm dot_xy_mm(2,:)];

%% Finish with z
z_mm = z_mm * ones(size(x_start_mm));