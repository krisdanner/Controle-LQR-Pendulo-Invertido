%==========================================================
% Controle por LQR de um pêndulo invertido com carrinho 
% e massa na ponta
%==========================================================
% Autor: Christian Danner Ramos de Carvalho
%==========================================================

close all;
clear;
clc;

% Define variáveis e parâmetros
syms x1 x2 x3 x4 u

M = 2;
m = 0.1;
b = 0.1;
l = 0.5;
g = 9.8;

% Desacopla as variáveis para obter as equações em função de dx e dtheta
U = [M+m m*l*cos(x3);
    cos(x3) l];

V = [u+m*l*x4^2*sin(x3);
    g*sin(x3)-b*x4];

W = inv(U)*V;

% Define a equação do pêndulo invertido com carrinho e força externa
f1 = x2;
f2 = W(1);
f3 = x4;
f4 = W(2);

% Linearize a equação em torno do ponto de equilíbrio para baixo (x3_eq=0)
% ou para cima (x3_eq=pi)
x1_eq = 0;
x2_eq = 0;
x3_eq = pi;
x4_eq = 0;
u_eq = 0;

A = jacobian([f1; f2; f3; f4], [x1; x2; x3; x4]);
A = subs(A, [x1, x2, x3, x4, u], [x1_eq, x2_eq, x3_eq, x4_eq, u_eq]);
A = double(A)

B = jacobian([f1; f2; f3; f4], u);
B = subs(B, [x1, x2, x3, x4, u], [x1_eq, x2_eq, x3_eq, x4_eq, u_eq]);
B = double(B)

C = [0 0 1 0];

D = 0;

sys = ss(A,B,C,D);

% Controlador LQR
Q = diag([1,1,100,1]);
R = 0.01;

K = lqr(A, B, Q, R);

% Sistema em malha fechada
    
Amf = A-B*K;
Bmf = B;
Cmf = C-D*K;
Dmf = D;

sysMF = ss(Amf,Bmf,Cmf,Dmf);

% Resposta do sistema em malha fechada 

t = 0:0.01:10;
u = zeros(size(t));
x0 = [0; 0; pi; 0];
[y,t,x] = lsim(sysMF,u,t,x0);
[y2,t2,x2] = lsim(sys,u,t,x0);

% Plot

figure;
plot(t,y,'linewidth',1);
grid;
hold on;
plot(t,y2,'linewidth',1);
legend('thetaMF','thetaMA')
xlabel('Tempo (s)');
ylabel('Ângulo do pêndulo (rad)');
title('Pêndulo invertido com controlador LQR');

% Plot da ação de controle
figure;
plot(t,-K*x','linewidth',1);
grid;
legend('u')
xlabel('Tempo (s)');
ylabel('Esforço de controle u');
title('Ação de controle');

% x0 = [10;0];
% 
% t = linspace(0,0.02,1000);
% [YMA tplotMA] = initial(planta,x0,t);
% [YMF tplotMF] = initial(plantaMF,x0,t);
% 
% figure;
% subplot(2,1,1);
% plot(t,YMA,'linewidth',2);
% xlabel('t(s)');
% ylabel('Vc(t)');
% legend('Malha Aberta');
% 
% subplot(2,1,2);
% plot(t,YMF,'linewidth',2);
% xlabel('t(s)');
% ylabel('Vc(t)');
% legend('Malha Fechada');


