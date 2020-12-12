%-------------------------------------------------------------------------%
%%Codigo Trabalho 2 feito por Mozart Fiorini Monteschio.
clear
clc
%-------------------------------------------------------------------------%
%%Declarando as Variáveis%%
J1 = 1;
J2 = 2;
K1 = 3;
K2 = 0.5;
b =0.25;
%-------------------------------------------------------------------------%
%% Parâmetros de simulação
to = 0;
tf = 30;
dt = 0.01;
t = [to:dt:tf];
N = numel(t);
%-------------------------------------------------------------------------%
%%Equações do Sistema%%
u = @(t) 10*sin(t);
u = u(t);
A = [-b/J1 -K1/J1 0 K1/J1; 1 0 0 0; 0 K1/J2 0 (-(K1+K2)/J2); 0 0 1 0];
B = [1/J1; 0; 0; 0];
C = [0 1 0 0; 0 0 0 1];
D = [0; 0];
%-------------------------------------------------------------------------%
%RK4
m0=[0; 0; 0; 0];
m = zeros(4, N);
art = zeros(2, N);
m(:, 1) = m0;
for k = 1:N-1
    art(:, k) = C * m(:, k) + D * entrada_u(dt * k); 
    m1 = A * m(:, k) + B * entrada_u(dt * k);
    m2 = A * (m(:, k) + dt/2 * m1) + B * entrada_u(dt * (k + 1/2));
    m3 = A * (m(:, k) + dt/2 * m2) + B * entrada_u(dt * (k + 1/2));
    m4 = A * (m(:, k) +  dt * m3) + B * entrada_u(dt * (k + 1));
    m(:, k + 1) = m(:, k) + dt * (m1/6 + m2/3 + m3/3 + m4/6);
end
%-------------------------------------------------------------------------%
plot(t, art)          %plot da saída
title('RK4')          %%Titulo do grafico 
ylabel('Saída')       %%Declarando nome do eixo Y
xlabel('Tempo')       %%Declarando nome do eixo X
legend('RK4 \theta_1', 'RK4 \theta_2')