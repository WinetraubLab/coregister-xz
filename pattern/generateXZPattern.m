function [...
    xStart_mm, xEnd_mm, ...
    yStart_mm, yEnd_mm, ...
    z_mm] = generateXZPattern()
% This function generates the instructions (lines) for XY lines
% The lines are from (x_start_mm(i), y_start_mm(i)) (to x_end_mm(i), y_end_mm(i)) 
% at depth z
% INPUT: verbose (default: false). If set to true, function will output
% line information as well as image 

%% Parameters
L_mm = 0.300;
D_mm = 0.20;
z0_mm = 40e-3; % [mm], how deep to draw pattern compared to the surface.
lensFOV_mm = 0.4;

xBuffer_mm = 25e-3; % This is to make sure that the diagonal lines don't
                    % mix with the up and down lines.

%% Init
xStart_mm=[];
xEnd_mm=[];
yStart_mm=[];
yEnd_mm =[];
z_mm = [];

txtPos_mm = []; % Position of the text for printing

%% Build

line_n = 0;

% Test different focal depths

generateZ(-0.825, -0.4:0.3:0.1, [0.07, 0.05, 0.03], false, L_mm*2);
%generateEdgePattern(-2*lensFOV_mm, -0.4:0.3:0.7, z0_mm + 0.075, false, L_mm);

generateZ(-0.575, -0.4, [0.08, 0.09, 0.10], true, L_mm*3);
%generateEdgePattern(-1.2*lensFOV_mm, -0.4:0.2:0.1, z0_mm + 0.050, true);


generateZ(-D_mm/2, -0.3:0.15:0.1, [z0_mm, z0_mm+0.025, z0_mm+0.025], true, L_mm);
xStart_mm(end)=[];xEnd_mm(end)=[];yStart_mm(end)=[];yEnd_mm(end)=[];z_mm(end)=[]; % Remove right most line
generateZ(+D_mm/2, -0.25:0.15:0.1, z0_mm + 0.000, true, L_mm);

generateZ(0.575, -0.3:0.15:0.1, [z0_mm + 0.000, z0_mm + 0.025, z0_mm + 0.025], true, L_mm);
%generateEdgePattern(+1.2*lensFOV_mm, -0.4:0.2:0.1, z0_mm + 0.050, true);
 
generateZ(0.825, -0.4:0.3:0.1, [z0_mm+0.000, z0_mm+0.000, z0_mm+0.000], false, L_mm*2);
%generateEdgePattern(+2*lensFOV_mm, -0.4:0.3:0.7, z0_mm + 0.075, false, L_mm);

%% Horizontal
xStart_mm = [xStart_mm(:)', -1.25];
xEnd_mm =   [xEnd_mm(:)',   +1.25];
yStart_mm = [yStart_mm(:)', -0.60];
yEnd_mm =   [yEnd_mm(:)',   -0.60];
z_mm = [z_mm(:)', z0_mm+00e-3];

%% Plot
figure(1)
for ii=1:length(xStart_mm)
    plot([xStart_mm(ii), xEnd_mm(ii)], [yStart_mm(ii), yEnd_mm(ii)]);
    if ii==1
        hold on;
    end
end
for ii=1:size(txtPos_mm,1)
    % text(txtPos_mm(ii,1),txtPos_mm(ii,2),sprintf('%d',ii-1),...
    %     'HorizontalAlignment','center',VerticalAlignment='middle');
    text(txtPos_mm(ii,1), txtPos_mm(ii,2), sprintf('%d', ii-1), ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize', 6);
end
plot(0.25*[-1 1 1 -1 -1], 0.25*[-1 -1 1 1 -1], 'k-', 'LineWidth', 2);
hold off
axis ij
axis equal

%% Generator functions
% This helper function generates a Z that is centered on x_0mm and starts
% at y0s_mm. If y0s_mm is an array - it will start at multiple positions
function generateZ(x0_mm, y0s_mm, z0_mm, isApplyLensYZone, L_mm)
    % Pattern looks like this:
    % |\  |
    % | \ | ^
    % |\ \| L
    % | \ | 
    % |  \|
    % < D >    
    
    % Up and down lines
    leftLine_x = x0_mm-D_mm/2;
    rightLine_x = x0_mm+D_mm/2;
    yMin = min(y0s_mm);
    yMax = max(y0s_mm) + L_mm;
    
    % Diagonal
    dx = D_mm;
    dy = L_mm;
    
    dStartX_mm = leftLine_x * ones(size(y0s_mm)) + xBuffer_mm;
    dStartY_mm = y0s_mm + dy/dx * xBuffer_mm;
    dEndX_mm   = rightLine_x * ones(size(y0s_mm)) - xBuffer_mm;
    dEndY_mm   = y0s_mm + L_mm - dy/dx * xBuffer_mm;
    
    % Assemble the line
    xStart1_mm = [leftLine_x, dStartX_mm(:)', rightLine_x];
    xEnd1_mm   = [leftLine_x, dEndX_mm(:)'  , rightLine_x];
    yStart1_mm = [yMin, dStartY_mm(:)', yMin];
    yEnd1_mm   = [yMax, dEndY_mm(:)'  , yMax];
    if isscalar(z0_mm)
        zLeftRail_mm = z0_mm;
        zRightRail_mm = z0_mm;
        zDiag_mm = ones(size(y0s_mm)) * z0_mm;
    elseif numel(z0_mm) == 3
        % [left rail, right rail, all diagonals]
        zLeftRail_mm = z0_mm(1);
        zDiag_mm = z0_mm(2);
        zRightRail_mm = ones(size(y0s_mm)) * z0_mm(3);
    elseif numel(z0_mm) == numel(y0s_mm)
        % Backward-compatible mode: one depth per diagonal, rails inherit edges.
        zDiag_mm = reshape(z0_mm, size(y0s_mm));
        zLeftRail_mm = zDiag_mm(1);
        zRightRail_mm = zDiag_mm(end);
    else
        error('generateZ:badDepthInput', ...
            ['z0_mm must be scalar, a 3-vector [left,right,diag], ' ...
             'or a vector matching the number of y0s_mm values.']);
    end

    
    % Normalize shapes for robust indexing regardless of row/column inputs.
    y0s_vec = y0s_mm(:)';
    zDiag_vec = zDiag_mm(:)';
    if numel(zDiag_vec) == 1 && numel(y0s_vec) > 1
        zDiag_vec = repmat(zDiag_vec, 1, numel(y0s_vec));
    end
    dStartX_vec = dStartX_mm(:)';
    dStartY_vec = dStartY_mm(:)';
    dEndX_vec = dEndX_mm(:)';
    dEndY_vec = dEndY_mm(:)';

    z1_mm = [zLeftRail_mm, zDiag_mm(:)', zRightRail_mm];

    % Apply zone
    if isApplyLensYZone
        ptStartIn = [xStart1_mm(:), yStart1_mm(:)]';
        ptEndIn = [xEnd1_mm(:), yEnd1_mm(:)]';
        [ptStart, ptEnd] = yOCTApplyEnableZone(ptStartIn,ptEndIn, ...
            @(x,y)(abs(y)<0.4/2));
        xStart1_mm = ptStart(1,:);
        yStart1_mm = ptStart(2,:);
        xEnd1_mm = ptEnd(1,:);
        yEnd1_mm = ptEnd(2,:);
    end

    % Print specs and add text labels for every generated line segment.
    % Left rail
    line_id = int32(line_n);
    z_left = zLeftRail_mm(1);
    fprintf("%d: {'line_type':'left_rail', 'L_mm':%.3f, 'D_mm':%.3f, 'x_offset_mm':%.3f, 'y_offset_mm':%.3f, 'z_mm':%.3f},\n",...
        line_id, L_mm, D_mm, x0_mm, yMin, z_left);
    txtPos_mm = [txtPos_mm; leftLine_x, mean([yMin yMax])];
    line_n = line_n+1;

    % Diagonals
    for i=1:numel(y0s_vec)
        line_id = int32(line_n);
        z_diag = zDiag_vec(i);
        fprintf("%d: {'line_type':'diag', 'L_mm':%.3f, 'D_mm':%.3f, 'x_offset_mm':%.3f, 'y_offset_mm':%.3f, 'z_mm':%.3f},\n",...
            line_id, L_mm, D_mm, x0_mm, y0s_vec(i), z_diag);
        txtPos_mm = [txtPos_mm; ...
            mean([dStartX_vec(i) dEndX_vec(i)]), mean([dStartY_vec(i) dEndY_vec(i)])];
        line_n = line_n+1;
    end

    % Right rail
    line_id = int32(line_n);
    z_right = zRightRail_mm(1);
    fprintf("%d: {'line_type':'right_rail', 'L_mm':%.3f, 'D_mm':%.3f, 'x_offset_mm':%.3f, 'y_offset_mm':%.3f, 'z_mm':%.3f},\n",...
        line_id, L_mm, D_mm, x0_mm, yMin, z_right);
    txtPos_mm = [txtPos_mm; rightLine_x, mean([yMin yMax])];
    line_n = line_n+1;

    % Append to lines
    xStart_mm = [xStart_mm(:)' xStart1_mm(:)'];
    xEnd_mm   = [xEnd_mm(:)'   xEnd1_mm(:)'];
    yStart_mm = [yStart_mm(:)' yStart1_mm(:)'];
    yEnd_mm   = [yEnd_mm(:)'   yEnd1_mm(:)'];
    z_mm   = [z_mm(:)' z1_mm(:)'];
end

% This helper function generates a || that helps locate which pattern
% are we on
function generateEdgePattern(x0,y0s,z0_mm, isApplyLensYZone, L_mm)
    
    yMin = min(y0s);
    if isApplyLensYZone
        yMin = max(yMin, -0.5/2);
    end
    yMax = max(y0s)+L_mm;
    if isApplyLensYZone
        yMax = min(yMax, +0.5/2);
    end
    yTotal_mm = yMax-yMin;
    yBin_mm = yTotal_mm/4;

    l = x0-D_mm/2-40e-3;
    r = x0+D_mm/2+40e-3;

    % Binarry pattern
    xStart_mm = [xStart_mm(:)' r            r            l];
    xEnd_mm   = [xEnd_mm(:)'   r            r            l];
    yStart_mm = [yStart_mm(:)' yMin+1*yBin_mm yMin+3*yBin_mm yMin+2*yBin_mm];
    yEnd_mm   = [yEnd_mm(:)'   yMin+2*yBin_mm yMin+4*yBin_mm yMin+4*yBin_mm];
    z_mm   = [z_mm(:)' ones(1,3)*z0_mm];
end
%% Main mission 
end
