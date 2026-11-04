function Ro = calculate_Ro(beta, mu, alpha, gamma, r)
    Ro = (beta * (mu + alpha)) / (mu * (mu + alpha) + gamma * (alpha * r + mu));
end