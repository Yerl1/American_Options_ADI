clc; clear all; close all;
M = 100;
N = 200;
T = 1; 
sigma = 0.3; 
r = 0.1; 
D0 = 0; 
X = 100; 

x_max = log(5); 
dx = x_max / M; 
x = linspace(0,x_max, M+1);

tau_exp = T * sigma^2 / 2; 
dtau = tau_exp / N; 


[S_f, P] = american_option(N,M);

S_f(end)
figure(1);
t =  linspace(0,1, N+1);
plot(t, S_f);

figure(2);
surf(t, x, P);

clc; clear all; close all;

% Параметры опциона
X = 100; 
r = 0.1; 
sigma = 0.3; 
T = 1; 
D0 = 0;

% Сетки для тестирования
M_values = [50, 100, 200]; % Количество шагов по x
N_values = [100, 200, 400]; % Количество шагов по времени

% Инициализация матриц для ошибок
RMSRE = zeros(length(M_values), length(N_values));

% Цикл по различным сеткам
for i = 1:length(M_values)
    for j = 1:length(N_values)
        M = M_values(i);
        N = N_values(j);
        
        % Численное решение
        [S_f_num, P_num] = american_option(N, M, X, r, sigma, T, D0);
        
        % Аналитическое решение Чжу (заглушка! Нужна ваша реализация)
        [S_f_analytical, P_analytical] = zhu_analytical_solution(N, M, X, r, sigma, T, D0);
        
        % Расчет относительной ошибки для S_f
        relative_error = abs(S_f_num - S_f_analytical) ./ S_f_analytical;
        RMSRE(i,j) = sqrt(mean(relative_error.^2));
    end
end

% Построение 3D-графика
figure;
[M_grid, N_grid] = meshgrid(M_values, N_values);
surf(M_grid, N_grid, RMSRE);
xlabel('(M)');
ylabel('(N)');
zlabel('RMSRE');
title(' RMSRE graph');