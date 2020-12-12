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
%Euler
x0=[0; 0; 0; 0];
x = zeros(4, N);
moz = zeros(2, N);
x(:, 1) = x0;
for k = 1:N-1
    moz(:, k) = C * x(:, k) + D * entrada_u(k * dt);
    x(:, k + 1) = (eye(4) + dt * A) * x(:, k)  + dt * B * entrada_u(k * dt); 
    
    
end
%-------------------------------------------------------------------------%
plot(t, moz)      %% %%Plotando Euler Explicito
title('Euler Explicito')   %%Titulo do grafico 
ylabel('Saída')       %%Declarando nome do eixo Y
xlabel('Tempo')       %%Declarando nome do eixo X
legend('Euler \theta_1','Euler \theta_2') %%Legenda para cada termo plotado
grid

