%-------------------------------------------------------------------------%
%%Codigo Trabalho 2 feito por Mozart Fiorini Monteschio. 18150371
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
A = [-b/J1 -K1/J1 0 K1/J1; 1 0 0 0; 0 K1/J2 0 (-(K1+K2)/J2); 0 0 1 0];
B = [1/J1; 0; 0; 0];
C = [0 1 0 0; 0 0 0 1];
D = [0; 0];
%-------------------------------------------------------------------------%
%Euler Explicito 
x0=[0; 0; 0; 0];
x = zeros(4, N);
moz = zeros(2, N);
x(:, 1) = x0;
for k = 1:N-1
    moz(:, k) = C * x(:, k) + D * entrada_u(k * dt);
    x(:, k +  1) = (eye(4) + dt * A) * x(:, k)  + dt * B * entrada_u(k * dt); 
end
%-------------------------------------------------------------------------%
% RK2 Metodo do Ponto Intermediario
a1 = 0;
a2 = 1;
p  = 1/2;
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
%-------------------------------------------------------------------------%
% RK2 Metodo do Heun
a11= 1/2;
a22 = 1/2;
p1  = 1;
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
%-------------------------------------------------------------------------%
%%RK2 Metodo do Ralston
a12 = 1/3;
a21 = 2/3;
p2  = 3/4;
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
%% Montando Espaços de Estado%%
sys = ss(A,B,C,D); 
sys = minreal(sys);
%-------------------------------------------------------------------------%
%%Comando Lsin%%
u = u(t);
y_lsim =lsim(sys,u,t);
na = size(t,2);
%-------------------------------------------------------------------------%
%%Plotando os Resultados em um unico grafico%%
plot(t,y_lsim,'r');       %%Plotando Lsin
hold on                   %%"Trava" o primeiro grafico para se plotar os outros por cima
plot(t, moz,'b')          %%Plotando Euler Explicito
plot(t, art,'g')          %%Plotando RK4 
plot(t,x(1,1:end-1),'k')  %%Plotando RK2 Metodo do Ponto Intermediario
plot(t,x2(1,1:end-1),'k') %%Plotando RK2 Metodo do Ponto Intermediario
plot(t,mot(1,1:end-1))    %%Plotando RK2 Metodo do Heun
plot(t,mot2(1,1:end-1))   %%Plotando RK2 Metodo do Heun
plot(t,mfm(1,1:end-1))    %%Plotando RK2 Metodo do Ralston
plot(t,mfm2(1,1:end-1))   %%Plotando RK2 Metodo do Ralston
title('Trabalho 2')       %%Titulo do grafico 
ylabel('Saída')           %%Declarando nome do eixo Y
xlabel('Tempo')           %%Declarando nome do eixo X
legend('Lsim \theta_1','Lsim \theta_2','Euler \theta_1','Euler \theta_2','RK4 \theta_1', 'RK4 \theta_2','RK2MPI \theta_1', 'RK2MPI \theta_2','RK2HEUN \theta_1', 'RK2HEUN \theta_2','RK2 RALSTON \theta_1', 'RK2 RALSTON \theta_2') %%Legenda para cada termo plotado
grid                      %%Coloca "Rede" atras do grafico,