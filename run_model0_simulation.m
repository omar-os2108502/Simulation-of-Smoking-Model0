function run_model0_simulation()
    % Clear workspace and active figures
    clc; clear; close all;

    % =========================================================================
    % PARAMETERS DIRECTLY FROM SHAROMI & GUMEL (2008), TABLE 2
    % =========================================================================
    mu    = 0.05;   % Natural mortality / turnover rate (per year)
    gamma = 0.30;   % Quitting rate (per year)
    alpha = 0.30;   % Relapse rate of temporary quitters (per year)
    r     = 0.50;   % Fraction quitting permanently
    N0    = 1660;   % Total constant population
    
    % Initial Conditions: [P(0), S(0), Qt(0), Qp(0)]
    y0    = [1600; 50; 10; 0]; 
    tspan = [0 150]; % 150 years
    
    % =========================================================================
    % SCENARIO 1: Eradication (Beta = 0.03 -> R0 < 1)
    % =========================================================================
    beta1 = 0.03;
    R0_1 = (beta1 * (mu + alpha)) / (mu * (mu + alpha) + gamma * (alpha * r + mu));
    [t1, y1] = ode45(@(t,y) odefun(t, y, beta1, mu, alpha, gamma, r, N0), tspan, y0);
    
    % Jacobian Matrix at Smoking-Free Equilibrium E0 = (N, 0, 0, 0)
    J_E0_scen1 = [ -mu, -beta1, 0, 0;
                     0,  beta1 - (mu + gamma), alpha, 0;
                     0,  gamma * (1 - r), -(mu + alpha), 0;
                     0,  gamma * r, 0, -mu ];
    eig1 = eig(J_E0_scen1);
    
    fprintf('=================================================================\n');
    fprintf('SCENARIO 1: Beta = %.2f | R0 = %.4f (< 1)\n', beta1, R0_1);
    fprintf('Eigenvalues at Smoking-Free Equilibrium E0:\n');
    disp(eig1);
    
    % =========================================================================
    % SCENARIO 2: Persistence / Endemic (Beta = 0.30 -> R0 > 1)
    % =========================================================================
    beta2 = 0.30;
    R0_2 = (beta2 * (mu + alpha)) / (mu * (mu + alpha) + gamma * (alpha * r + mu));
    [t2, y2] = ode45(@(t,y) odefun(t, y, beta2, mu, alpha, gamma, r, N0), tspan, y0);
    
    % Jacobian Matrix at Smoking-Free Equilibrium E0 = (N, 0, 0, 0)
    J_E0_scen2 = [ -mu, -beta2, 0, 0;
                     0,  beta2 - (mu + gamma), alpha, 0;
                     0,  gamma * (1 - r), -(mu + alpha), 0;
                     0,  gamma * r, 0, -mu ];
    eig2 = eig(J_E0_scen2);
    
    % Numerical Steady-State Endemic Equilibrium E*
    P_star  = y2(end, 1);
    S_star  = y2(end, 2);
    Qt_star = y2(end, 3);
    Qp_star = y2(end, 4);
    
    % Jacobian Matrix at Endemic Equilibrium E*
    J_Estar = [ -(beta2 * S_star / N0) - mu, -(beta2 * P_star / N0), 0, 0;
                 (beta2 * S_star / N0),       (beta2 * P_star / N0) - (mu + gamma), alpha, 0;
                 0,                          gamma * (1 - r), -(mu + alpha), 0;
                 0,                          gamma * r, 0, -mu ];
    eig_Estar = eig(J_Estar);
    
    fprintf('=================================================================\n');
    fprintf('SCENARIO 2: Beta = %.2f | R0 = %.4f (> 1)\n', beta2, R0_2);
    fprintf('Eigenvalues at Smoking-Free Equilibrium E0:\n');
    disp(eig2);
    fprintf('Endemic Equilibrium E* = (P=%.1f, S=%.1f, Qt=%.1f, Qp=%.1f)\n', P_star, S_star, Qt_star, Qp_star);
    fprintf('Eigenvalues at Endemic Equilibrium E*:\n');
    disp(eig_Estar);
    fprintf('=================================================================\n');

    % =========================================================================
    % PLOTTING TRAJECTORIES
    % =========================================================================
    figure('Color', 'k', 'Position', [100, 100, 1100, 480]);
    
    % Subplot 1: Eradication
    subplot(1, 2, 1);
    plot(t1, y1(:,1), 'b-',  'LineWidth', 1.8); hold on;
    plot(t1, y1(:,2), 'r--', 'LineWidth', 1.8);
    plot(t1, y1(:,3), 'm-.', 'LineWidth', 1.8);
    plot(t1, y1(:,4), 'g:',  'LineWidth', 1.8);
    title(sprintf('Eradication: R_0 = %.4f < 1', R0_1), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Time (Years)', 'FontSize', 11);
    ylabel('Population', 'FontSize', 11);
    legend('Potential Smokers (P)', 'Active Smokers (S)', 'Temp Quitters (Q_t)', 'Perm Quitters (Q_p)', 'Location', 'east');
    grid on;
    
    % Subplot 2: Persistence
    subplot(1, 2, 2);
    plot(t2, y2(:,1), 'b-',  'LineWidth', 1.8); hold on;
    plot(t2, y2(:,2), 'r--', 'LineWidth', 1.8);
    plot(t2, y2(:,3), 'm-.', 'LineWidth', 1.8);
    plot(t2, y2(:,4), 'g:',  'LineWidth', 1.8);
    title(sprintf('Persistence: R_0 = %.4f > 1', R0_2), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Time (Years)', 'FontSize', 11);
    ylabel('Population', 'FontSize', 11);
    legend('Potential Smokers (P)', 'Active Smokers (S)', 'Temp Quitters (Q_t)', 'Perm Quitters (Q_p)', 'Location', 'east');
    grid on;
end

% ODE System
function dydt = odefun(~, y, beta, mu, alpha, gamma, r, N)
    P  = y(1);
    S  = y(2);
    Qt = y(3);
    Qp = y(4);
    
    dP  = mu * N - (beta * P * S / N) - mu * P;
    dS  = (beta * P * S / N) + alpha * Qt - (mu + gamma) * S;
    dQt = gamma * (1 - r) * S - (mu + alpha) * Qt;
    dQp = gamma * r * S - mu * Qp;
    
    dydt = [dP; dS; dQt; dQp];
end