classdef VerticalLine
    % Defines a vertical (rail) line segment.
    %
    % Parameters
    % ----------
    % x       : X position of the line [mm]
    % yStart  : Y start position [mm]
    % yEnd    : Y end position [mm]
    % depth   : Z depth of the line [mm]

    properties
        x      (1,1) double
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
        function obj = VerticalLine(x, yStart, yEnd, depth)
            obj.x      = x;
            obj.yStart = yStart;
            obj.yEnd   = yEnd;
            obj.depth  = depth;
        end

        function val = get.xStart_mm(obj)
            val = obj.x;
        end

        function val = get.xEnd_mm(obj)
            val = obj.x;
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