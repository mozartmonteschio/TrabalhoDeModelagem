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
a12 = 1/3;
a21 = 2/3;
p2  = 3/4;
%-------------------------------------------------------------------------%
% RK2 Metodo do Ralston
mfm = zeros(2,length(t)); %pre alocação t1,w1
mfm(:,1)=[0; 0]; % Condic¸oes iniciais t1,w1
mfm2 = zeros(2,length(t)); % pre alocação t2,w2
mfm2(:,1)=[0; 0]; % Condic¸oes iniciais t2,w2
for k=1:length(t)
Vin(k) =  u(t(k)); % Armazenamento do vetor de entrada
m11(k) =  mfm(2,k); %teta1
m12(k) =  mfm2(2,k); %teta2
m13(k) =  -b/J1*mfm(2,k)-K1/J1*mfm(1,k)+K1/J1*mfm2(1,k)+u(t(k)/J1); %w1   
m14(k) =  K1/J2*mfm(1,k) - K2/J2*mfm2(1,k) - K1/J2*mfm2(1,k); %w2
m21(k) =  mfm(2,k)+dt*p2*m13(k); %W1 + dt*p*m13
m22(k) =  mfm2(2,k)+dt*p2*m14(k);; %teta 2
m23(k) =  (-b/J1*mfm(2,k)+dt*p2*m13(k))-(K1/J1*mfm(1,k)+dt*p2*m11(k))+(K1/J1*mfm2(1,k)+dt*p2*m12(k))+u(t(k)+p2*dt/J1);
m24(k) =  (K1/J2*mfm(1,k)+dt*p2*m11(k)) - (K2/J2*mfm2(1,k)+dt*p2*m12(k)) - (K1/J2*mfm2(1,k)+dt*p2*m12(k));
mfm(1,k+1) = mfm(1,k)+dt*(a12*m11(k)+a21*m21(k)); %teta1
mfm(2,k+1) = mfm(2,k)+dt*(a12*m13(k)+a21*m23(k)); %w1
mfm2(1,k+1)= mfm2(1,k)+dt*(a12*m12(k)+a21*m22(k)); %teta2
mfm2(2,k+1)= mfm2(2,k)+dt*(a12*m14(k)+a21*m24(k)); %w2
end
plot(t,mfm(1,1:end-1))
hold on
plot(t,mfm2(1,1:end-1))
title('RK2 Metodo do Ralston')   %%Titulo do grafico 
ylabel('Saída')       %%Declarando nome do eixo Y
xlabel('Tempo')       %%Declarando nome do eixo X
legend('\theta_1','\theta_2') %%Legenda para cada termo plotado
grid                  %%Coloca "Rede" atras do grafico,