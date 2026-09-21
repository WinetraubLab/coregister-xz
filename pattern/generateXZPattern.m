function [xStart_mm, xEnd_mm, yStart_mm, yEnd_mm, z_mm] = generateXZPattern()
% Generates the instructions (lines) for XY lines.
% Constructors:
%   VerticalLine(x, yStart, yEnd, depth)
%   DiagonalLine(xStart, xEnd, yStart, yEnd, depth)

%% Build line objects
allLines = {};

allLines{end+1} = VerticalLine(-0.9, 0.5, -0.5, 65e-3);
allLines{end+1} = DiagonalLine(-0.88, -0.77, -0.5, 0.5, 65e-3);
allLines{end+1} = VerticalLine(-0.75, 0.5, -0.5, 40e-3);

% center left
allLines{end+1} = VerticalLine(-0.525, -0.375, 0.375, 40e-3); % longer
allLines{end+1} = DiagonalLine(-0.5, -0.39, 0.02, -0.17, 65e-3);
allLines{end+1} = DiagonalLine(-0.5, -0.39, 0.17, -0.02, 65e-3);
allLines{end+1} = VerticalLine(-0.375, 0.17, -0.17, 65e-3);
allLines{end+1} = DiagonalLine(-0.35, -0.25, 0.02, -0.17, 40e-3);
allLines{end+1} = DiagonalLine(-0.35, -0.25, 0.17, -0.02, 40e-3);
allLines{end+1} = VerticalLine(-0.225, 0.17, -0.17, 65e-3);

% center 
allLines{end+1} = VerticalLine(0.025, 0.17, -0.17, 65e-3);
allLines{end+1} = DiagonalLine(0.027, 0.15, -0.17, 0.02, 40e-3);
allLines{end+1} = DiagonalLine(0.027, 0.15, -0.02, 0.17, 40e-3);
allLines{end+1} = VerticalLine(0.175, 0.17, -0.17, 65e-3);
allLines{end+1} = DiagonalLine(0.19, 0.3, -0.17, 0.02, 65e-3);
allLines{end+1} = DiagonalLine(0.19, 0.3, -0.02, 0.17, 65e-3);
allLines{end+1} = VerticalLine(0.325, -0.17, 0.17, 40e-3); 

% center right
allLines{end+1} = VerticalLine(0.575, 0.17, -0.17, 65e-3);
allLines{end+1} = DiagonalLine(0.59, 0.7, -0.17, 0.02, 40e-3);
allLines{end+1} = DiagonalLine(0.59, 0.7, -0.02, 0.17, 40e-3);
allLines{end+1} = VerticalLine(0.725, 0.17, -0.17, 65e-3);
allLines{end+1} = DiagonalLine(0.75, 0.85, -0.17, 0.02, 65e-3);
allLines{end+1} = DiagonalLine(0.75, 0.85, -0.02, 0.17, 65e-3);
allLines{end+1} = VerticalLine(0.875, -0.375, 0.375, 40e-3); % longer

allLines{end+1} = VerticalLine(1.1, 0.5, -0.5, 40e-3);
allLines{end+1} = DiagonalLine(1.12, 1.23, 0.5, -0.5, 65e-3);
allLines{end+1} = VerticalLine(1.25, 0.5, -0.5, 65e-3);

% straight line across the top
allLines{end+1} = DiagonalLine(-1, 1.7, 0.5, 0.5, 65e-3);

%% Convert to arrays
[xStart_mm, xEnd_mm, yStart_mm, yEnd_mm, z_mm] = linesToArrays(allLines);

%% Preview
plotPattern(allLines, xStart_mm, xEnd_mm, yStart_mm, yEnd_mm);

end

function [xStart_mm, xEnd_mm, yStart_mm, yEnd_mm, z_mm] = linesToArrays(lines)
% Unpack line objects into parallel coordinate arrays.

n = numel(lines);

xStart_mm = zeros(1, n);
xEnd_mm   = zeros(1, n);
yStart_mm = zeros(1, n);
yEnd_mm   = zeros(1, n);
z_mm      = zeros(1, n);

for ii = 1:n
    xStart_mm(ii) = lines{ii}.xStart_mm;
    xEnd_mm(ii)   = lines{ii}.xEnd_mm;
    yStart_mm(ii) = lines{ii}.yStart_mm;
    yEnd_mm(ii)   = lines{ii}.yEnd_mm;
    z_mm(ii)      = lines{ii}.z_mm;
end

end

function plotPattern(allLines, xStart_mm, xEnd_mm, yStart_mm, yEnd_mm)
%Draw lines with index labels and the OCT scan area outlined.

% ROI box overlaid on the pattern
oct_xRange = [-0.525 0.875];
oct_yRange = [-0.175 0.175];

figure(1); clf;
hold on;

for ii = 1:numel(allLines)
    plot([xStart_mm(ii), xEnd_mm(ii)], [yStart_mm(ii), yEnd_mm(ii)]);

    x_mid = (xStart_mm(ii) + xEnd_mm(ii)) / 2;
    y_mid = (yStart_mm(ii) + yEnd_mm(ii)) / 2;
    text(x_mid, y_mid, sprintf('%d', ii - 1), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 8, ...
        'FontWeight', 'bold');
end

plot(oct_xRange([1 2 2 1 1]), ...
     oct_yRange([1 1 2 2 1]), ...
     'k-', 'LineWidth', 2);

hold off;
axis ij;
axis equal;

end
