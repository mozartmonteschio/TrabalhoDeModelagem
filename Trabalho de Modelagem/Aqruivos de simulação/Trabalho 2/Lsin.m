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
%% Montando Espaços de Estado%%
sys = ss(A,B,C,D); 
sys = minreal(sys);
%-------------------------------------------------------------------------%
%%Comando Lsin%%
y_lsim =lsim(sys,u,t);
na = size(t,2);
%-------------------------------------------------------------------------%
plot(t,y_lsim);      
title('Lsim')   
ylabel('Saída')      
xlabel('Tempo')      
legend('Lsim \theta_1','Lsim \theta_2') %%Legenda para cada termo plotado
grid                  %%Coloca "Rede" atras do grafico,