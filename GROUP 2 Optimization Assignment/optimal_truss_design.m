function optimal_truss_design()
    P = 20000;     % Load in N
    sigma_0 = 150e6; % Allowable stress in Pa
    rho = 7850;    % Density kg/m^3
    l = 1.0;       % Base length m

    % Objective: Minimize weight
    obj_fun = @(x) rho * l * (2*sqrt(2)*x(1) + x(2));

    nonlin_con = @(x) truss_constraints(x, P, sigma_0);

    x0 = [0.001, 0.001]; % Initial areas in m^2
    lb = [1e-5, 1e-5];
    ub = [0.01, 0.01];

    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    [x_opt, fval] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, nonlin_con, options);

    fprintf('\n--- Optimal Three-Bar Truss Solution ---\n');
    fprintf('Optimal Area of Outer Bars (A1): %.6f m^2\n', x_opt(1));
    fprintf('Optimal Area of Inner Bar (A2):  %.6f m^2\n', x_opt(2));
    fprintf('Minimum Structural Weight:       %.4f kg\n', fval);

    % --- Plotting ---
    figure;
    [A1, A2] = meshgrid(linspace(0.0001, 0.003, 100), linspace(0.0001, 0.003, 100));
    Weight = rho * l .* (2*sqrt(2)*A1 + A2);
    
    % Stress constraints evaluation
    s1 = P * (sqrt(2)*A1 + A2) ./ (A1.^2 + sqrt(2)*A1.*A2);
    s2 = P * (sqrt(2)*A1 + 2*A2) ./ (sqrt(2)*A1.^2 + 2*A1.*A2);
    s3 = P * (A2 - sqrt(2)*A1) ./ (A1.^2 + sqrt(2)*A1.*A2);
    
    Weight(s1 > sigma_0 | s2 > sigma_0 | s3 > sigma_0) = NaN;

    contourf(A1, A2, Weight, 30, 'LineColor', 'none');
    colorbar;
    hold on;
    plot(x_opt(1), x_opt(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'r');
    xlabel('Area $A_1$ (m^2)', 'Interpreter', 'latex');
    ylabel('Area $A_2$ (m^2)', 'Interpreter', 'latex');
    title('Feasible Region & Optimum for Three-Bar Truss Design');
    grid on;
end

function [c, ceq] = truss_constraints(x, P, sigma_0)
    A1 = x(1);
    A2 = x(2);
    
    c(1) = (P * (sqrt(2)*A1 + A2)) / (A1^2 + sqrt(2)*A1*A2) - sigma_0;
    c(2) = (P * (sqrt(2)*A1 + 2*A2)) / (sqrt(2)*A1^2 + 2*A1*A2) - sigma_0;
    c(3) = (P * (A2 - sqrt(2)*A1)) / (A1^2 + sqrt(2)*A1*A2) - sigma_0;
    ceq = [];
end