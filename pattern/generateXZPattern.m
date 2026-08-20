function [xStart_mm, xEnd_mm, yStart_mm, yEnd_mm, z_mm] = generateXZPattern()
% Generates the instructions (lines) for XY lines.
% Lines go from (xStart_mm(i), yStart_mm(i)) to (xEnd_mm(i), yEnd_mm(i)) at depth z.

%% Build line objects
allLines = {};

% VerticalLine(x, yStart, yEnd, depth)
% DiagonalLine(xStart, xEnd, yStart, yEnd, depth)

allLines{end+1} = VerticalLine(-1, 0.5, -0.5, 65e-3);
allLines{end+1} = VerticalLine(-0.92, 0.5, -0.5, 40e-3);

allLines{end+1} = VerticalLine(-0.7, -0.375, 0.375, 40e-3);
allLines{end+1} = DiagonalLine(-0.67, -0.53, 0.37, -0.37, 65e-3);
allLines{end+1} = VerticalLine(-0.5, -0.375, 0.375, 65e-3);

% in FOV:
allLines{end+1} = VerticalLine(-0.25, 0.25, -0.25, 40e-3);
allLines{end+1} = DiagonalLine(-0.22, -0.03, -0.25, -0.05, 65e-3);
allLines{end+1} = DiagonalLine(-0.22, -0.03, -0.1, 0.1, 65e-3);
allLines{end+1} = DiagonalLine(-0.22, -0.03, 0.05, 0.25, 65e-3);

allLines{end+1} = VerticalLine(0.0, 0.25, -0.25, 65e-3);
allLines{end+1} = DiagonalLine(0.03, 0.22, -0.2, 0.0, 40e-3);
allLines{end+1} = DiagonalLine(0.03, 0.22, -0.05, 0.15, 40e-3);
allLines{end+1} = DiagonalLine(0.03, 0.22, 0.1, 0.3, 40e-3);

allLines{end+1} = VerticalLine(0.25, 0.25, -0.25, 65e-3);

% out of FOV
allLines{end+1} = VerticalLine(0.45, 0.375, -0.375, 65e-3);
allLines{end+1} = DiagonalLine(0.48, 0.62, 0.375, -0.13, 65e-3);
allLines{end+1} = DiagonalLine(0.48, 0.62, 0.13, -0.375, 65e-3);
allLines{end+1} = VerticalLine(0.65, 0.375, -0.375, 40e-3);

allLines{end+1} = VerticalLine(0.9, 0.5, -0.5, 65e-3);
allLines{end+1} = VerticalLine(1.1, 0.5, -0.5, 40e-3);
allLines{end+1} = DiagonalLine(0.93, 1.07, 0.5, -0.5, 40e-3);



%% Convert to arrays
[xStart_mm, xEnd_mm, yStart_mm, yEnd_mm, z_mm] = linesToArrays(allLines);

%% Plot
figure(1); clf;
for ii = 1:numel(allLines)
    plot([xStart_mm(ii), xEnd_mm(ii)], [yStart_mm(ii), yEnd_mm(ii)]);
    if ii == 1; hold on; end
end
fov_half = 0.25;
plot(fov_half*[-1 1 1 -1 -1], fov_half*[-1 -1 1 1 -1], 'k-', 'LineWidth', 2);
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
