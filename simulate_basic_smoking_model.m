function simulate_basic_smoking_model()

%% 1. Parameter Set
mu = 1/20;    % Natural death rate
alpha = 0.3;  % Rate temporary quitters revert to somking
gamma = 0.3;  % Quitting rate
r = 0.5;      % Fraction who quit permanently
% Intial population setup: [P, S, Qt, Qp]
P0  = 1600;
S0  = 50;
Qt0 = 10;
Qp0 = 0;
N0  = P0 + S0 + Qt0 + Qp0;
y0  = [P0; S0; Qt0; Qp0];
    
tspan = [0, 60]; % Simulation time in years
%% 2. Scenario 1: Rs < 1 (e.g., beta = 0.03)
    beta1 = 0.03;
    Rs1 = calculate_Rs(beta1, mu, alpha, gamma, r);
    [t1, y1] = ode45(@(t, y) basic_smoking_ode(t, y, beta1, mu, alpha, gamma, r, N0), tspan, y0);
    
    %% 3. Scenario 2: Rs > 1 (e.g., beta = 0.3)
    beta2 = 0.3;
    Rs2 = calculate_Rs(beta2, mu, alpha, gamma, r);
    [t2, y2] = ode45(@(t, y) basic_smoking_ode(t, y, beta2, mu, alpha, gamma, r, N0), tspan, y0);
    
    %% 4. Plotting Results
    figure('Name', 'Basic Smoking Model Simulation', 'Color', 'k');
    
    % Plot 1: Rs < 1
    subplot(1, 2, 1);
    plot(t1, y1(:,1), 'b', 'LineWidth', 1.5); hold on;
    plot(t1, y1(:,2), 'r', 'LineWidth', 1.5);
    plot(t1, y1(:,3), 'm', 'LineWidth', 1.5);
    plot(t1, y1(:,4), 'g', 'LineWidth', 1.5);
    grid on;
    title(sprintf('R_s = %.4f < 1 (Smoking Eradicated)', Rs1));
    xlabel('Time (Years)'); ylabel('Population');
    legend('Potential (P)', 'Smokers (S)', 'Temp Quit (Q_t)', 'Perm Quit (Q_p)');
    
    % Plot 2: Rs > 1
    subplot(1, 2, 2);
    plot(t2, y2(:,1), 'b', 'LineWidth', 1.5); hold on;
    plot(t2, y2(:,2), 'r', 'LineWidth', 1.5);
    plot(t2, y2(:,3), 'm', 'LineWidth', 1.5);
    plot(t2, y2(:,4), 'g', 'LineWidth', 1.5);
    grid on;
    title(sprintf('R_s = %.4f > 1 (Smoking Persists)', Rs2));
    xlabel('Time (Years)'); ylabel('Population');
    legend('Potential (P)', 'Smokers (S)', 'Temp Quit (Q_t)', 'Perm Quit (Q_p)');
end