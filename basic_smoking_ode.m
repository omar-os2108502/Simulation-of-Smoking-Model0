function dydt = basic_smoking_ode(~, y, beta, mu, alpha, gamma, r, N)
    P  = y(1);
    S  = y(2);
    Qt = y(3);
    Qp = y(4);
    
    dP  = mu*N - (beta * P * S / N) - mu * P;               % Eq (1)
    dS  = (beta * P * S / N) + alpha * Qt - (mu + gamma)*S; % Eq (2)
    dQt = gamma * (1 - r) * S - (mu + alpha) * Qt;          % Eq (3)
    dQp = gamma * r * S - mu * Qp;                          % Eq (4)
    
    dydt = [dP; dS; dQt; dQp];
end