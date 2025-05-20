function S_f = american_option(N, M)
T = 1;
sigma = 0.3;
r = 0.1;
D0 = 0;
X = 100;

x_max = log(5);
dx = x_max / M;
x = 0:dx:x_max;

tau_exp = T * sigma^2 / 2;
dtau = tau_exp / N;
tau = 0:dtau:tau_exp;

gamma = 2 * r / sigma^2;
D = 2 * D0 / sigma^2;

P = zeros(M+1, N+1);
S_f = zeros(1, N+1);
S_f(1) = 1;

for n = 1:N
    P_curr = P(:, n);
    S_f_curr = S_f(n);

    % --- Predictor Step ---
    alpha = 1 + (gamma/2) * dx^2;
    beta = 1 + dx + ((D + 1)/2) * dx^2;
    S_f_hat = (alpha - P_curr(2)) / beta; % Predict S_f_hat

    % --- Corrector Step ---
    % Compute S_f_ratio for coefficients
    S_f_ratio = (S_f_hat - S_f_curr) / (S_f_hat + S_f_curr);

    % Coefficients for matrix A and B
    a = (dtau/(2*dx^2)) - (gamma - D -1)*(dtau/(4*dx)) - (1/(2*dx)) * S_f_ratio;
    c = (dtau/(2*dx^2)) + (gamma - D -1)*(dtau/(4*dx)) + (1/(2*dx)) * S_f_ratio;
    b = 1 + dtau/dx^2 + 0.5*gamma*dtau;
    b_prime = 1 - dtau/dx^2 - 0.5*gamma*dtau;

    % Construct A and B matrices
    main_diag = b * ones(M-1, 1);
    lower_diag = -a * ones(M-2, 1);
    upper_diag = -c * ones(M-2, 1);

    lower_diag = [0; lower_diag];
    upper_diag = [upper_diag; 0];

    A = spdiags([lower_diag, main_diag, upper_diag], [-1, 0, 1], M-1, M-1);

    main_diag_B = b_prime * ones(M-1, 1);
    lower_diag_B = a * ones(M-2, 1);
    upper_diag_B = c * ones(M-2, 1);

    lower_diag_B = [0; lower_diag_B];
    upper_diag_B = [upper_diag_B; 0];
    B = spdiags([lower_diag_B, main_diag_B, upper_diag_B], [-1, 0, 1], M-1, M-1);

    % Boundary terms: P₀ⁿ and P₀^{n+1} (use S_f_curr and S_f_hat)
    boundary_term = a * ( (1 - S_f_curr) + (1 - S_f_hat) );
    RHS = B * P_curr(2:M) + boundary_term * [1; zeros(M-2, 1)];

    % Solve linear system
    P_new_inner = A \ RHS;
    P(2:M, n+1) = P_new_inner;

    % --- Update S_f(n+1) using corrected P₁ ---
    S_f(n+1) = (alpha - P(2, n+1)) / beta;

    % --- Corrected Boundary Condition: Update P₀^{n+1} using S_f(n+1) ---
    P(1, n+1) = 1 - S_f(n+1); % This was the critical fix
    P(M+1, n+1) = 0;
end

% Convert normalized S_f back to actual price
S_f_real = S_f * X;
end
