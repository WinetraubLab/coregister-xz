function [xStart_mm, xEnd_mm, yStart_mm, yEnd_mm, z_mm] = generateXZPattern()
% Generates the instructions (lines) for XY lines.
% Lines go from (xStart_mm(i), yStart_mm(i)) to (xEnd_mm(i), yEnd_mm(i)) at depth z.

%% Build line objects
allLines = {};

% VerticalLine(x, yStart, yEnd, depth)
% DiagonalLine(xStart, xEnd, yStart, yEnd, depth)

allLines{end+1} = VerticalLine(-0.9, 0.5, -0.5, 65e-3);
allLines{end+1} = DiagonalLine(-0.88, -0.77, -0.5, 0.5, 65e-3);
allLines{end+1} = VerticalLine(-0.75, 0.5, -0.5, 40e-3);

% center left:
allLines{end+1} = VerticalLine(-0.525, -0.375, 0.375, 40e-3); % longer
allLines{end+1} = DiagonalLine(-0.5, -0.36, 0.02, -0.17, 65e-3);
allLines{end+1} = DiagonalLine(-0.5, -0.36, 0.17, -0.02, 65e-3);
allLines{end+1} = VerticalLine(-0.33, 0.17, -0.17, 65e-3);
allLines{end+1} = DiagonalLine(-0.30, -0.17, 0.02, -0.17, 40e-3);
allLines{end+1} = DiagonalLine(-0.30, -0.17, 0.17, -0.02, 40e-3);
allLines{end+1} = VerticalLine(-0.15, 0.17, -0.17, 65e-3);

% center right:
allLines{end+1} = VerticalLine(0.15, 0.17, -0.17, 65e-3);
allLines{end+1} = DiagonalLine(0.17, 0.33, -0.17, 0.02, 40e-3);
allLines{end+1} = DiagonalLine(0.17, 0.33, -0.02, 0.17, 40e-3);
allLines{end+1} = VerticalLine(0.33, 0.17, -0.17, 65e-3);
allLines{end+1} = DiagonalLine(0.36, 0.5, -0.17, 0.02, 65e-3);
allLines{end+1} = DiagonalLine(0.36, 0.5, -0.02, 0.17, 65e-3);
allLines{end+1} = VerticalLine(0.525, -0.375, 0.375, 40e-3); % longer

allLines{end+1} = VerticalLine(0.75, 0.5, -0.5, 40e-3);
allLines{end+1} = DiagonalLine(0.88, 0.77, 0.5, -0.5, 65e-3);
allLines{end+1} = VerticalLine(0.9, 0.5, -0.5, 65e-3);

allLines{end+1} = DiagonalLine(-1, 1, 0.5, 0.5, 65e-3);

%% Convert to arrays
[xStart_mm, xEnd_mm, yStart_mm, yEnd_mm, z_mm] = linesToArrays(allLines);

%% Plot
figure(1); clf;
hold on;
for ii = 1:numel(allLines)
    plot([xStart_mm(ii), xEnd_mm(ii)], [yStart_mm(ii), yEnd_mm(ii)]);

    % Midpoint of each line for the label
    x_mid = (xStart_mm(ii) + xEnd_mm(ii)) / 2;
    y_mid = (yStart_mm(ii) + yEnd_mm(ii)) / 2;
    text(x_mid, y_mid, sprintf('%d', ii-1), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 8, ...
        'FontWeight', 'bold');
end
% fov_half = 0.25;
fov_half_x = 0.525;
fov_half_y = 0.176;

plot(fov_half_x*[-1 1 1 -1 -1], ...
     fov_half_y*[-1 -1 1 1 -1], ...
     'k-', 'LineWidth', 2);

hold off
axis ij
axis equal
end

function [xStart_mm, xEnd_mm, yStart_mm, yEnd_mm, z_mm] = linesToArrays(lines)

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
