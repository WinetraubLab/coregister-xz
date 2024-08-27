% This script will generate a 3D volume of how the XZ pattern should appear
% in the tissue

% Load the pattern
[x_start_mm, x_end_mm, y_start_mm, y_end_mm, z_mm] = ...
    generateXZPattern();

%% Script Inputs

% Define the bounding box to generate 3D simulation for
x_grid_mm = -1.00:2e-3:1.00;
y_grid_mm = -1.00:2e-3:1.00;
z_grid_mm =  0.00:2e-3:0.10;

% Physical parameters 
NA = 0.35; % Match NA to observed photobleach pattern. For 40x use 0.35 (though lens NA is 0.8)
lambda_mm = 900e-9*1e3; % Wavelength in m
n = 1.4; % Medium index of refraction
photobleach_intensity = 40/(diff(y_grid_mm(1:2))/1e-3); % Can be any number >0

% Plot OCT volume on top?
oct_scan_mm = [-0.25 0.25]; % OCT x-y scan size
% If you don't want to plot oct_scan_mm, un-comment:
% oct_scan_mm = nan;

% Simulation output
output_tiff_file = 'out_xz.tiff';

%% Configurable Parameters
% Gausian base waist
w0_mm = 1/pi*lambda_mm*n/NA;
zR_mm = pi*w0_mm^2/lambda_mm;

%% Create a gread
[xx_mm, yy_mm] = meshgrid(x_grid_mm,y_grid_mm);
pixel_size_mm = diff(x_grid_mm(1:2));
output_data = zeros(length(z_grid_mm),length(x_grid_mm),length(y_grid_mm),'uint8');

% Loop for each plane in z, start at the bottom to match coordinate system
% orientation: https://docs.google.com/presentation/d/1tOod76WvhvOuByo-K81YB4b3QjRq-6A5j2ztS_ANSNo/edit#slide=id.g2812f1bd1d1_0_52
isFirstLoop=true;
for zi=length(z_grid_mm):-1:1
    z=z_grid_mm(zi);
    c_all = ones(size(xx_mm));

    % Loop over all lines
    for lineI=1:length(x_start_mm)

        if y_start_mm(lineI) == y_end_mm(lineI)
            % Horizontal line  
            yI = abs(yy_mm-y_start_mm(lineI)) < pixel_size_mm;
            xI = xx_mm >= min(x_start_mm(lineI),x_end_mm(lineI)) & ...
                 xx_mm <= max(x_start_mm(lineI),x_end_mm(lineI));

        elseif x_start_mm(lineI) == x_end_mm(lineI)
            % Vertical line
            xI = abs(xx_mm-x_start_mm(lineI))<pixel_size_mm;
            yI = yy_mm >= min(y_start_mm(lineI),y_end_mm(lineI)) & ...
                 yy_mm <= max(y_start_mm(lineI),y_end_mm(lineI));
        end


        % Gausian waist at the depth we are
        wz_mm = w0_mm*sqrt(1+( ...
            (z_mm(lineI)-z) /zR_mm)^2);

        % Add the line
        c=createPhotobleachArea(yI&xI,photobleach_intensity,wz_mm/pixel_size_mm);
        c_all = c_all .* c;
    end

    % Add the OCT scan
    if all(~isnan(oct_scan_mm))
        c_all = addOCTScanRectangle(c_all, xx_mm, yy_mm, oct_scan_mm,pixel_size_mm);
    end

    output_data(zi,:,:) = c_all;
    
    % Present to user
    figure(27);
    imagesc(x_grid_mm,y_grid_mm,c_all)
    axis equal
    caxis([0,1])
    colormap gray
    pause(0.1);
end

%% Save
yOCT2Tif(output_data,output_tiff_file);

%% Finally, print to screen a description of how this pattern was created
describePattern(x_start_mm, x_end_mm, y_start_mm, y_end_mm, z_mm)

function c_all = addOCTScanRectangle(c_all, xx_mm, yy_mm, oct_scan_mm,pixel_size_mm)
r = zeros(size(c_all));

% Left and right borders
r(yy_mm >= oct_scan_mm(1) & yy_mm <= oct_scan_mm(2) & abs(xx_mm-oct_scan_mm(1)) < pixel_size_mm) = 1;
r(yy_mm >= oct_scan_mm(1) & yy_mm <= oct_scan_mm(2) & abs(xx_mm-oct_scan_mm(2)) < pixel_size_mm) = 1;

% Top and bottom borders 
r(xx_mm >= oct_scan_mm(1) & xx_mm <= oct_scan_mm(2) & abs(yy_mm-oct_scan_mm(1)) < pixel_size_mm) = 1;
r(xx_mm >= oct_scan_mm(1) & xx_mm <= oct_scan_mm(2) & abs(yy_mm-oct_scan_mm(2)) < pixel_size_mm) = 1;

c_all(r==1) = 0;

end

function describePattern(x_start_mm, x_end_mm, y_start_mm, y_end_mm, z_mm)
z = unique(z_mm);

fprintf('Draw Lines from (x,y) to (x,y) at the following depths. Units: mm\n');
for i=1:length(z)
    fprintf('Depth: %.3f\n',z(i))
    ii = find(z_mm == z(i));
    for j=ii
        fprintf('   (%+.3f,%+.3f) -> (%+.3f,%+.3f)\n',...
            x_start_mm(j),y_start_mm(j),...
            x_end_mm(j),y_end_mm(j)...
            );
    end
end
end

function c=createPhotobleachArea(whereToPhotobleach,photobleach_intensity,w)
% whereToPhotobleach - bolean mask = 1 if photobleach center / seed, 0
%   otherwise.
% w - gausian waist size in pixels

% Create canvas 
c = ones(size(whereToPhotobleach))*(photobleach_intensity+1); 

% Create the photobleach area
c(whereToPhotobleach) = 0; 

% Gausian filt
c = imgaussfilt(c, w/sqrt(2))-photobleach_intensity;

% Clean uo
c(c<0)=0;

end
