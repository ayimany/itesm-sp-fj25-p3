%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Definitions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function points = draw_parametric_wire(parametric_function, current_function, t_range, step)
    if (t_range(1) > t_range(2))
        error("t-range min cannot be greater than t_range max")
    end

    point_count = round((t_range(2) - t_range(1)) / step);
    points = zeros(point_count, 4);

    i = 0;
    for t = t_range(1):step:t_range(2)
        i = i + 1;
        point = subs(parametric_function, t);
        current = subs(current_function, t);
        points(i, 1) = point(1);
        points(i, 2) = point(2);
        points(i, 3) = point(3);
        points(i, 4) = current;
    end
end

function [magnetic_field_x, magnetic_field_y, magnetic_field_z] = calculate_magnetic_field_influence(points, x_grid, y_grid, z_grid)
    mu_0 = 4 * pi * 1E-7;

    magnetic_field_x = zeros(size(x_grid));
    magnetic_field_y = zeros(size(y_grid));
    magnetic_field_z = zeros(size(z_grid));
    
    for i = 1:size(points, 1)-1
        point_start = points(i, 1:3);
        point_end = points(i+1, 1:3);
        current = points(i, 4);
        
        segment_delta = point_end - point_start;
        segment_midpoint = (point_start + point_end) / 2;
        
        x_segment_displacement = x_grid - segment_midpoint(1);
        y_segment_displacement = y_grid - segment_midpoint(2);
        z_segment_displacement = z_grid - segment_midpoint(3);
        
        displacement_magnitude = sqrt(x_segment_displacement.^2 + y_segment_displacement.^2 + z_segment_displacement.^2);
        displacement_magnitude(displacement_magnitude == 0) = eps;
        
        % Vector section of the Biot-Savart law
        cross_x = (segment_delta(2) * z_segment_displacement - segment_delta(3) * y_segment_displacement);
        cross_y = (segment_delta(3) * x_segment_displacement - segment_delta(1) * z_segment_displacement);
        cross_z = (segment_delta(1) * y_segment_displacement - segment_delta(2) * x_segment_displacement);
        
        % Non vector section of the Biot-Savart law
        coeff = mu_0 * current / (4 * pi) ./ (displacement_magnitude.^2);
        
        magnetic_field_x = magnetic_field_x + coeff .* cross_x;
        magnetic_field_y = magnetic_field_y + coeff .* cross_y;
        magnetic_field_z = magnetic_field_z + coeff .* cross_z;
    end
end

function curves = add_parametric_curve(fx, fy, fz, fi, range, curves)
    curves = [curves; fx, fy, fz, fi, range];
end

%%%%%%%%%%%%%%%%%
%%% CONSTANTS %%%
%%%%%%%%%%%%%%%%%

syms t;

simulation_limit = 2.5;
simulation_step = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SPACES AND RANGES %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

[x_axis, y_axis, z_axis] = meshgrid(-simulation_limit : simulation_step : simulation_limit);

%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER DEFINITIONS %%%
%%%%%%%%%%%%%%%%%%%%%%%%

parametric_curves = [];

parametric_curves = add_parametric_curve(t    , 0    , 0, 1, [0, 1], parametric_curves);
parametric_curves = add_parametric_curve(1    , t - 1, 0, 1, [1, 2], parametric_curves);
parametric_curves = add_parametric_curve(3 - t, 1    , 0, 1, [2, 3], parametric_curves);
parametric_curves = add_parametric_curve(0    , 4 - t, 0, 1, [3, 4], parametric_curves);

%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETRIC WIRES %%%
%%%%%%%%%%%%%%%%%%%%%%%%

parametric_wires = [];

for i = 1 : height(parametric_curves)
    p_obj = parametric_curves(i, :);
    curve = [p_obj(1), p_obj(2), p_obj(3)];
    current = p_obj(4);
    r_lo = double(p_obj(5));
    r_hi = double(p_obj(6));

    parametric_wires = cat(3, parametric_wires, draw_parametric_wire(curve, current, [r_lo, r_hi], simulation_step));
end

%%%%%%%%%%%%%%%%%%%%%%
%%% MAGNETIC FIELD %%%
%%%%%%%%%%%%%%%%%%%%%%

magnetic_field_x = zeros(size(x_axis));
magnetic_field_y = zeros(size(y_axis));
magnetic_field_z = zeros(size(z_axis));

for i = 1 : size(parametric_wires, 3)
    wire = parametric_wires(:, :, i);
    [bx, by, bz] = calculate_magnetic_field_influence(wire, x_axis, y_axis, z_axis);
    magnetic_field_x = magnetic_field_x + bx;
    magnetic_field_y = magnetic_field_y + by;
    magnetic_field_z = magnetic_field_z + bz;
end

%%%%%%%%%%%%%%%%
%%% PLOTTING %%%
%%%%%%%%%%%%%%%%

figure;
hold on;
view(3);

title('Magnetic Field');
xlabel('x'); ylabel('y'); zlabel('z');

grid on;
axis equal;

quiver3(x_axis, y_axis, z_axis, magnetic_field_x, magnetic_field_y, magnetic_field_z);

for i = 1 : size(parametric_wires, 3)
    wire = parametric_wires(:, :, i);
    plot3(wire(:, 1), wire(:, 2), wire(:, 3), '-r');
end

