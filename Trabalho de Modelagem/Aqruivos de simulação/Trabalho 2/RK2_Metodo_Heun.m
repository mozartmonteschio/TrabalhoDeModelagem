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
a11= 1/2;
a22 = 1/2;
p1  = 1;
%-------------------------------------------------------------------------%
% RK2 Metodo do Heun
mot = zeros(2,length(t)); %pre alocação t1,w1
mot(:,1)=[0; 0]; % Condic¸oes iniciais t1,w1
mot2 = zeros(2,length(t)); % pre alocação t2,w2
mot2(:,1)=[0; 0]; % Condic¸oes iniciais t2,w2
for k=1:length(t)
Vin(k) =  u(t(k)); % Armazenamento do vetor de entrada
m11(k) =  mot(2,k); %teta1
m12(k) =  mot2(2,k); %teta2
m13(k) =  -b/J1*mot(2,k)-K1/J1*mot(1,k)+K1/J1*mot2(1,k)+u(t(k)/J1); %w1   
m14(k) =  K1/J2*mot(1,k) - K2/J2*mot2(1,k) - K1/J2*mot2(1,k); %w2
m21(k) =  mot(2,k)+dt*p1*m13(k); %W1 + dt*p*m13
m22(k) =  mot2(2,k)+dt*p1*m14(k);; %teta 2
m23(k) =  (-b/J1*mot(2,k)+dt*p1*m13(k))-(K1/J1*mot(1,k)+dt*p1*m11(k))+(K1/J1*mot2(1,k)+dt*p1*m12(k))+u(t(k)+p1*dt/J1);
m24(k) =  (K1/J2*mot(1,k)+dt*p1*m11(k)) - (K2/J2*mot2(1,k)+dt*p1*m12(k)) - (K1/J2*mot2(1,k)+dt*p1*m12(k));
mot(1,k+1) = mot(1,k)+dt*(a11*m11(k)+a22*m21(k)); %teta1
mot(2,k+1) = mot(2,k)+dt*(a11*m13(k)+a22*m23(k)); %w1
mot2(1,k+1)= mot2(1,k)+dt*(a11*m12(k)+a22*m22(k)); %teta2
mot2(2,k+1)= mot2(2,k)+dt*(a11*m14(k)+a22*m24(k)); %w2
end
plot(t,mot(1,1:end-1))
hold on
plot(t,mot2(1,1:end-1))
title('RK2 Metodo do Heun')   %%Titulo do grafico 
ylabel('Saída')       %%Declarando nome do eixo Y
xlabel('Tempo')       %%Declarando nome do eixo X
legend('\theta_1','\theta_2') %%Legenda para cada termo plotado
grid                  %%Coloca "Rede" atra