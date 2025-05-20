clc; clear all; close all;
M = 100;
N = 100;
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

[S_f, P] = american_option(N,M);
S_f = S_f * 100;


disp(S_f(end))

figure(1);
plot(x, S_f, 'b');


figure(2);
surf(tau, x, P);