classdef DiagonalLine
    % Defines a diagonal line segment.
    %
    % Parameters
    % ----------
    % xStart  : X start position [mm]
    % xEnd    : X end position [mm]
    % yStart  : Y start position [mm]
    % yEnd    : Y end position [mm]
    % depth   : Z depth of the line [mm]

    properties
        xStart (1,1) double
        xEnd   (1,1) double
        yStart (1,1) double
        yEnd   (1,1) double
        depth  (1,1) double
    end

    properties (Dependent)
        xStart_mm
        xEnd_mm
        yStart_mm
        yEnd_mm
        z_mm
    end

    methods
        function obj = DiagonalLine(xStart, xEnd, yStart, yEnd, depth)
            obj.xStart = xStart;
            obj.xEnd   = xEnd;
            obj.yStart = yStart;
            obj.yEnd   = yEnd;
            obj.depth  = depth;
        end

        function val = get.xStart_mm(obj)
            val = obj.xStart;
        end

        function val = get.xEnd_mm(obj)
            val = obj.xEnd;
        end

        function val = get.yStart_mm(obj)
            val = obj.yStart;
        end

        function val = get.yEnd_mm(obj)
            val = obj.yEnd;
        end

        function val = get.z_mm(obj)
            val = obj.depth;
        end
    end
end