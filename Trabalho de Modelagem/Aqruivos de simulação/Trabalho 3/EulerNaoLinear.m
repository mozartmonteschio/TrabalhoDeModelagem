%-------------------------------------------------------------------------%
%%Codigo Trabalho 2 feito por Mozart Fiorini Monteschio. 18150371
clear
clc
%-------------------------------------------------------------------------%
%%Equa��es do Sistema e como condicoes iniciais%%
A =[0 1;-20 -2];
B = [0;0];
C = [1 0];
D = [0];
%-------------------------------------------------------------------------%
%% Par�metros de simula��o
to = 0;
tf = 12;
dt = 0.01;
t = [to:dt:tf];
na=size(t,2);
xo=[(3*pi)/4;0]; 
%-------------------------------------------------------------------------%
%Simulacao euler n�o linear
x(:,1)=xo;
for k=1:na
m1 = x(1,k);
m2 = x(2,k);
x(1, k+1) = m2*dt+m1;
x(2, k+1) = m2+dt*(-2*m2-20*sin(m1));
x1(k+1)=x(1, k+1);
end
%-------------------------------------------------------------------------%
%%Plotando os Resultados
plot(t,x1(:,1:end-1),'g')
title('Euler Simul��o 1')
ylabel('Saida')
xlabel('Tempo')
grid
