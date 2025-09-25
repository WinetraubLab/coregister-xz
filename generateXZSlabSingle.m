function [x_start_mm, x_end_mm, y_start_mm, y_end_mm, z_mm] = ...
    generateXZSlabSingle(L_sq_mm, pitch_mm, z_list_mm, centerX_mm, centerY_mm, orientation, includeBorders)
% Genera UN solo cuadrado (slab) de lado L_sq_mm con líneas paralelas
% y devuelve vectores de líneas + z por línea.
%
% Entradas:
%   L_sq_mm        - lado del cuadrado (mm), p.ej. 0.50
%   pitch_mm       - espaciado entre líneas (mm), p.ej. 0.02 (= 20 µm)
%   z_list_mm      - profundidad(es) (mm). Escalar o vector [z1 z2 ...]
%   centerX_mm     - centro X del cuadrado
%   centerY_mm     - centro Y del cuadrado
%   orientation    - 'v' (vertical) o 'h' (horizontal)
%   includeBorders - true/false para incluir los dos bordes del cuadrado
%
% Salidas:
%   x_start_mm, x_end_mm, y_start_mm, y_end_mm, z_mm

    if nargin < 7, includeBorders = true; end
    if nargin < 6 || isempty(orientation), orientation = 'v'; end

    half   = L_sq_mm/2;
    x_left = centerX_mm - half;  x_right = centerX_mm + half;
    y_top  = centerY_mm - half;  y_bot   = centerY_mm + half;

    x_start_mm = []; x_end_mm = []; y_start_mm = []; y_end_mm = []; z_mm = [];

    switch orientation
        case 'v' % líneas verticales (X fijo, Y recorre)
            xs = [];
            if includeBorders
                xs = [xs, x_left, x_right];
            end
            n_interior = max(0, floor(L_sq_mm/pitch_mm) - 1);
            if n_interior > 0
                xs_int = linspace(x_left + pitch_mm, x_right - pitch_mm, n_interior);
                xs = [xs, xs_int];
            end
            x_start_mm = xs;
            x_end_mm   = xs;
            y_start_mm = y_top  * ones(size(xs));
            y_end_mm   = y_bot  * ones(size(xs));

        case 'h' % líneas horizontales (Y fijo, X recorre)
            ys = [];
            if includeBorders
                ys = [ys, y_top, y_bot];
            end
            n_interior = max(0, floor(L_sq_mm/pitch_mm) - 1);
            if n_interior > 0
                ys_int = linspace(y_top + pitch_mm, y_bot - pitch_mm, n_interior);
                ys = [ys, ys_int];
            end
            x_start_mm = x_left  * ones(size(ys));
            x_end_mm   = x_right * ones(size(ys));
            y_start_mm = ys;
            y_end_mm   = ys;

        otherwise
            error('orientation debe ser ''v'' o ''h''.');
    end

    % Asignación de profundidades
    nLines = numel(x_start_mm);
    if isscalar(z_list_mm)
        z_mm = z_list_mm * ones(1, nLines);
    else
        z_mm = zeros(1, nLines);
        for i = 1:nLines
            z_mm(i) = z_list_mm(1 + mod(i-1, numel(z_list_mm)));
        end
    end
end
