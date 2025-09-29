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
D_mm = 0.150;
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

generateZ(-2*lensFOV_mm, -0.4:0.3:0.1, z0_mm + 0.000, false, L_mm*2);
%generateEdgePattern(-2*lensFOV_mm, -0.4:0.3:0.7, z0_mm + 0.075, false, L_mm);

generateZ(-1*lensFOV_mm, -0.3:0.15:0.1, z0_mm + 0.025, true, L_mm);
%generateEdgePattern(-1.2*lensFOV_mm, -0.4:0.2:0.1, z0_mm + 0.050, true);

generateZ(-D_mm/2, -0.3:0.15:0.1, z0_mm + 0.000, true, L_mm);
xStart_mm(end)=[];xEnd_mm(end)=[];yStart_mm(end)=[];yEnd_mm(end)=[];z_mm(end)=[]; % Remove right most line
generateZ(+D_mm/2, -0.25:0.15:0.1, z0_mm + 0.000, true, L_mm);

generateZ(+1*lensFOV_mm, -0.3:0.15:0.1, z0_mm + 0.025, true, L_mm);
%generateEdgePattern(+1.2*lensFOV_mm, -0.4:0.2:0.1, z0_mm + 0.050, true);

generateZ(+2*lensFOV_mm, -0.4:0.3:0.1, z0_mm + 0.000, false, L_mm*2);
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
    text(txtPos_mm(ii,1),txtPos_mm(ii,2),sprintf('%d',ii-1),...
        'HorizontalAlignment','center',VerticalAlignment='middle');
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
    z1_mm = ones(size(xStart1_mm))*z0_mm;

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

    % Print the line spec & add a text for it
    for i=1:length(y0s_mm)
        fprintf("%d: {'L_mm':%.3f, 'D_mm':%.3f, 'x_offset_mm':%.3f, 'y_offset_mm':%.3f, 'z_mm':%.3f},\n",...
            line_n, L_mm,D_mm,x0_mm,y0s_mm(i),z0_mm)
        
        txtPos_mm = [txtPos_mm; ...
            mean([leftLine_x rightLine_x]), mean([dStartY_mm(i) dEndY_mm(i)])];

        line_n = line_n+1;
    end
    
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
