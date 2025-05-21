function [S_f,P] = american_option(N, M)
T = 1;
sigma = 0.3;
X = 100;
r = 0.1;
D0 = 0;

gamma = 2*r/sigma^2;
D = 2*D0/sigma^2;

x_max = log(5);
dx = x_max/M;
x = linspace(0, x_max, M+1);


tau_exp = T * sigma^2 / 2;
dtau = tau_exp/N; 
tau = linspace(0, tau_exp, N+1);

P = zeros(M+1, N+1);
S_f = zeros(N+1,1);
S_f(1) = 1;

alpha = 1 + (gamma/2)*dx^2;
beta = 1 + dx + (D + 1)*dx^2/2;

for i = 1:N

    P_curr = P(:, i);

    % Predictor 
    S_f_hat = (alpha - P_curr(2)) / beta;

    % Corrector 
    a_coeff = dtau/(2*dx^2) - (gamma - D - 1)*dtau/(4*dx) - (S_f_hat - S_f(i))/(2*dx*(S_f_hat + S_f(i)));
    b_coeff = 1 + dtau/dx^2 + gamma*dtau/2;
    c_coeff = dtau/(2*dx^2) + (gamma - D - 1)*dtau/(4*dx) + (S_f_hat - S_f(i))/(2*dx*(S_f_hat + S_f(i)));

    main_diag = b_coeff * ones(M-1, 1);
    lower_diag = -a_coeff * ones(M-1, 1);
    upper_diag = -c_coeff * ones(M-1, 1);

    A = spdiags([lower_diag, main_diag, upper_diag], -1:1, M-1, M-1);
    B_main = (1 - dtau/dx^2 - gamma*dtau/2) * ones(M-1, 1);
    B = spdiags([a_coeff*ones(M-1,1), B_main, c_coeff*ones(M-1,1)], -1:1, M-1, M-1);

    RHS = B * P_curr(2:M) + [a_coeff*(1 - S_f(i) + (1 - S_f_hat)); zeros(M-3,1); c_coeff*P_curr(M)];
    P_next_inner = A \ RHS;

    % Update P and S_f
    P(2:M, i+1) = P_next_inner;
    S_f(i+1) = (alpha - P(2, i+1)) / beta;
    P(1, i+1) = 1 - S_f(i+1); 
    P(M+1, i+1) = 0;          

end

S_f = S_f * X;
end

% 
% 
% function [S_f,P] = american_option(N, M)
% T = 1;
% sigma = 0.3;
% X = 100;
% r = 0.1;
% D0 = 0;
% 
% gamma = 2*r/sigma^2;
% D = 2*D0/sigma^2;
% 
% x_max = log(5); 
% dx = x_max/M;
% x = linspace(0, x_max, M+1);
% 
% 
% tau_exp = T * sigma^2 / 2; 
% dtau = tau_exp/N; 
% tau = linspace(0, tau_exp, N+1);
% 
% P = zeros(M+1, N+1); 
% S_f = zeros(N+1,1);
% S_f(1) = 1;  
% 
% alpha = 1 + (gamma/2)*dx^2;
% beta = 1 + dx + (D + 1)*dx^2/2;
% 
% for n = 1:N
%     P_curr = P(:,n);
%     % Predictor
%     C = (P_curr(3) - 2 * P_curr(2) + P_curr(1))/dx^2;
%     termD = (gamma - D - 1) * (P_curr(3) - P_curr(1))/(2*dx);
%     E = gamma * P_curr(2);
%     % H = (P(2,n+1) - P_curr(2))/dtau;
%     F = (P_curr(3) - P_curr(1))/(2*dx*S_f(n));
%     S_f_hat = (alpha - P_curr(2) + F * S_f(n) - (C + termD - E)* dtau)/(F + beta);
%     % S_f_hat = (H - C - termD + E) * dtau/ F + S_f(n);
% 
%     % Corrector
%     a_coeff = dtau/(2*dx^2) - (gamma - D - 1)*dtau/(4*dx) - (S_f_hat - S_f(n))/(2*dx*(S_f_hat + S_f(n)));
%     b_coeff = 1 + dtau/dx^2 + gamma*dtau/2;
%     c_coeff = dtau/(2*dx^2) + (gamma - D - 1)*dtau/(4*dx) + (S_f_hat - S_f(n))/(2*dx*(S_f_hat + S_f(n)));
% 
%     main_diag = b_coeff * ones(M-1, 1);
%     lower_diag = -a_coeff * [ones(M-2, 1); 0];
%     upper_diag = -c_coeff * [0; ones(M-2, 1)];
% 
%     A = spdiags([lower_diag, main_diag, upper_diag], -1:1, M-1, M-1);
%     B_main = (1 - dtau/dx^2 - gamma*dtau/2) * ones(M-1, 1);
%     B = spdiags([a_coeff*ones(M-1,1), B_main, c_coeff*ones(M-1,1)], -1:1, M-1, M-1);
% 
%     RHS = B * P_curr(2:M) + [a_coeff * (P_curr(1) + (1 - S_f_hat)); zeros(M-3,1); c_coeff * P_curr(end)];
%     P_next_inner = A \ RHS;
% 
%     P(2:M, n+1) = P_next_inner;
%     S_f(n+1) = (alpha - P(2, n+1)) / beta;
%     P(1, n+1) = 1 - S_f(n+1); 
% end
% S_f = S_f * X;
% end