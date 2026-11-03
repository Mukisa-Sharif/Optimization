function optimal_beam_design()
    w = 15000;         % Distributed load in N/m
    L = 4.0;           % Span in m
    E = 200e9;         % Young's modulus in Pa
    sigma_allow = 166.7e6; % Allowable bending stress in Pa
    delta_allow = L / 360; % Max allowable deflection in m

    % Objective function: Minimize Area (b * h)
    obj_fun = @(x) x(1) * x(2);

    % Constraints
    nonlin_con = @(x) beam_constraints(x, w, L, E, sigma_allow, delta_allow);

    % Initial guess [b, h]
    x0 = [0.1, 0.3];

    % Bounds: 0.02 <= b <= 0.5, 0.05 <= h <= 1.0
    lb = [0.02, 0.05];
    ub = [0.5, 1.0];

    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    [x_opt, fval] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, nonlin_con, options);

    fprintf('\n--- Optimal Rectangular Beam Solution ---\n');
    fprintf('Optimal Width (b): %.4f m\n', x_opt(1));
    fprintf('Optimal Height (h): %.4f m\n', x_opt(2));
    fprintf('Minimum Cross-Sectional Area: %.4f m^2\n', fval);

    % --- Plotting ---
    figure;
    [B, H] = meshgrid(linspace(0.05, 0.3, 100), linspace(0.1, 0.6, 100));
    Area = B .* H;
    
    M_max = (w * L^2) / 8;
    sigma_max = M_max ./ ((B .* H.^2) ./ 6);
    I_beam = (B .* H.^3) ./ 12;
    delta_max = (5 * w * L^4) ./ (384 * E * I_beam);
    
    Area(sigma_max > sigma_allow | delta_max > delta_allow | (H./B < 1.5) | (H./B > 3.5)) = NaN;

    contourf(B, H, Area, 25, 'LineColor', 'none');
    colorbar;
    hold on;
    plot(x_opt(1), x_opt(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'r');
    xlabel('Beam Width $b$ (m)', 'Interpreter', 'latex');
    ylabel('Beam Height $h$ (m)', 'Interpreter', 'latex');
    title('Feasible Region & Optimum for Rectangular Beam Design');
    grid on;
end

function [c, ceq] = beam_constraints(x, w, L, E, sigma_allow, delta_allow)
    b = x(1);
    h = x(2);
    
    M_max = (w * L^2) / 8;
    Z = (b * h^2) / 6;
    I = (b * h^3) / 12;
    
    sigma = M_max / Z;
    delta = (5 * w * L^4) / (384 * E * I);
    
    c(1) = sigma - sigma_allow;          % Stress constraint
    c(2) = delta - delta_allow;          % Deflection constraint
    c(3) = 1.5 * b - h;                  % Lower aspect ratio (h >= 1.5b)
    c(4) = h - 3.5 * b;                  % Upper aspect ratio (h <= 3.5b)
    ceq = [];
end