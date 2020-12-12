%-------------------------------------------------------------------------%
%%Codigo Trabalho 2 feito por Mozart Fiorini Monteschio. 18150371
clear
clc
%-------------------------------------------------------------------------%
%%Equa��es do Sistema e como condicoes iniciais%%
A =[0 1;-20 -2];;
B = [0;4];
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
%Simula�ao euler linear
u = ones(1,na);
x2(:,1)=xo;
for k =1:na 
  x2(:,k+1) = (A*dt + eye(2))*x2(:,k);
end
%-------------------------------------------------------------------------%
%%Plotando os Resultados
plot(t,x2(1,1:end-1),'b')
title('Euler Simul��o 2')
xlabel('Tempo')
ylabel('Saida')
grid