function optimal_shaft_design()
% Parameters
T = 12000;         % Torque in N*m
L = 1.5;           % Length in m
G = 79e9;          % Shear modulus in Pa
tau_allow = 130e6; % Allowable shear stress in Pa
theta_allow = 1.5 * (pi/180); % Allowable angle of twist in rad

%  % Objective function: Minimize Area
    obj_fun = @(x) pi * (x(1)^2 - x(2)^2);

    % Non-linear constraints function
    nonlin_con = @(x) shaft_constraints(x, T, L, G, tau_allow, theta_allow);

    % Initial guess [ro, ri]
    x0 = [0.08, 0.05];

    % Bounds: lb <= x <= ub
    lb = [0.01, 0.001];
    ub = [0.20, 0.18];

    % Optimization using fmincon
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    [x_opt, fval] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, nonlin_con, options);

    fprintf('\n--- Optimal Hollow Shaft Solution ---\n');
    fprintf('Optimal Outer Radius (ro): %.4f m\n', x_opt(1));
    fprintf('Optimal Inner Radius (ri): %.4f m\n', x_opt(2));
    fprintf('Minimum Cross-Sectional Area: %.4f m^2\n', fval);

    % --- Plotting Feasible Region & Optimum ---
    figure;
    [R_o, R_i] = meshgrid(linspace(0.04, 0.12, 100), linspace(0.01, 0.10, 100));
    Area = pi .* (R_o.^2 - R_i.^2);
    
    % Evaluate constraints across grid for contour plotting
    J_grid = (pi/2) .* (R_o.^4 - R_i.^4);
    stress_viol = (T .* R_o ./ J_grid) - tau_allow;
    twist_viol = (T * L ./ (G .* J_grid)) - theta_allow;
    geom_viol = R_i - 0.9 .* R_o;
    
    % Mask infeasible points
    Area(stress_viol > 0 | twist_viol > 0 | geom_viol > 0) = NaN;

    contourf(R_o, R_i, Area, 30, 'LineColor', 'none');
    colorbar;
    hold on;
    plot(x_opt(1), x_opt(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'r');
    xlabel('Outer Radius $r_o$ (m)', 'Interpreter', 'latex');
    ylabel('Inner Radius $r_i$ (m)', 'Interpreter', 'latex');
    title('Feasible Design Space & Optimum for Hollow Shaft');
    grid on;
end

function [c, ceq] = shaft_constraints(x, T, L, G, tau_allow, theta_allow)
    ro = x(1);
    ri = x(2);
    J = (pi/2) * (ro^4 - ri^4);
    
    % c <= 0 form
    c(1) = (T * ro / J) - tau_allow;             % Stress constraint
    c(2) = (T * L / (G * J)) - theta_allow;      % Twist constraint
    c(3) = ri - 0.9 * ro;                        % Thickness ratio constraint
    ceq = [];
end
