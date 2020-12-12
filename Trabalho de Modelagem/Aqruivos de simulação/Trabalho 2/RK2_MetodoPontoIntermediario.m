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
b = 0.25;
%-------------------------------------------------------------------------%
%% Parâmetros de simulação
to = 0;
tf = 30;
dt = 0.01;
dt = 0.01;
t = [to:dt:tf];
N = numel(t);
%-------------------------------------------------------------------------%
%%Equações do Sistema%%
u  = @(t) 10*sin(t);
a1 = 0;
a2 = 1;
p  = 1/2;
%-------------------------------------------------------------------------%
% RK2 Metodo do Ponto Intermediario
x = zeros(2,length(t)); %pre alocação t1,w1
x(:,1)=[0; 0]; % Condic¸oes iniciais t1,w1
x2 = zeros(2,length(t)); % pre alocação t2,w2
x2(:,1)=[0; 0]; % Condic¸oes iniciais t2,w2
for k=1:length(t)
Vin(k) =  u(t(k)); % Armazenamento do vetor de entrada
m11(k) =  x(2,k); %teta1
m12(k) =  x2(2,k); %teta2
m13(k) =  -b/J1*x(2,k)-K1/J1*x(1,k)+K1/J1*x2(1,k)+u(t(k)/J1); %w1   
m14(k) =  K1/J2*x(1,k) - K2/J2*x2(1,k) - K1/J2*x2(1,k); %w2
m21(k) =  x(2,k)+dt*p*m13(k); %W1 + dt*p*m13
m22(k) =  x2(2,k)+dt*p*m14(k);; %teta 2
m23(k) =  (-b/J1*x(2,k)+dt*p*m13(k))-(K1/J1*x(1,k)+dt*p*m11(k))+(K1/J1*x2(1,k)+dt*p*m12(k))+u(t(k)+p*dt/J1);
m24(k) =  (K1/J2*x(1,k)+dt*p*m11(k)) - (K2/J2*x2(1,k)+dt*p*m12(k)) - (K1/J2*x2(1,k)+dt*p*m12(k));
x(1,k+1) = x(1,k)+dt*(a1*m11(k)+a2*m21(k)); %teta1
x(2,k+1) = x(2,k)+dt*(a1*m13(k)+a2*m23(k)); %w1
x2(1,k+1)= x2(1,k)+dt*(a1*m12(k)+a2*m22(k)); %teta2
x2(2,k+1)= x2(2,k)+dt*(a1*m14(k)+a2*m24(k)); %w2
end
plot(t,x(1,1:end-1))
hold on
plot(t,x2(1,1:end-1))
title('RK2 Metodo do Ponto Intermediario')   %%Titulo do grafico 
ylabel('Saída')       %%Declarando nome do eixo Y
xlabel('Tempo')       %%Declarando nome do eixo X
legend('\theta_1','\theta_2') %%Legenda para cada termo plotado
grid                  %%Coloca "Rede" atras do grafico,