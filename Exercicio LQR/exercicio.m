close all; clc; clear;
%% Modelo massa-mola-amortecedor triplo
m1=0.5; %kg
m2=0.5; %kg
m3=0.5; %kg
k1=1;  %N/cm
k2=1;  %N/cm
b1=0.02; %Ns/cm
b2=0.01; %Ns/cm
A=[   0         0        0      1         0        0   ; ...
      0         0        0      0         1        0   ; ...
      0         0        0      0         0        1   ; ...
   -k1/m1    +k1/m1      0   -b1/m1    +b1/m1      0   ; ...
   +k1/m2 -(k1+k2)/m2 +k2/m2 +b1/m2 -(b1+b2)/m2 +b2/m2 ; ...
      0      +k2/m3   -k2/m3    0      +b2/m3   -b2/m3   ...
];
B=[   0     0   ;...
      0     0   ;...
      0     0   ;...
    +1/m1   0   ;...
      0     0   ;...
      0   -1/m3  ...
 ];
C=eye(6);
Cx=[eye(3) zeros(3)];
Cv=[zeros(3) eye(3)];
D=zeros(6,2);
%% Matrizes de ponderação do LQR
% Regra de Bryson
Qii=[20;20;20;10;10;10].^-2;
Rii=[20 20].^-2;
Q=diag(Qii);
R=diag(Rii);
x0=[-5 5 5 1 0 -1]';
%% Controlabilioade
Cc = ctrb(A,B);
%% Calculo do ganho LQR
%% LQR de Horizonte Finito 
T=3; 
PT=zeros(6); % Control effort is zero at t = T
% PT = eye(6); % Other designed scenario
% Solução da Equação Diferencial de Riccati (numérica)
opt = odeset('AbsTol',1.0e-8,'RelTol',1.0e-8);
[tvP,vP] = ode45(@(t,x) riccatiDE(t,x,A,B,Q,R), [T,0], reshape(PT,6*6,1),opt);
% Solução do problema de condiçãop inicial com o controle LQR
[t,x] = ode45(@(t,x) MMAmodel(t,x,A,B,R,tvP,vP), [0,T], x0,opt);
%% LQR de Horizonto Infinito
% Pc=care(A,B,Q,R); % Solução da Equação Algébrica e Riccati
% Klqr = R\B'*P; % Ganho LQR de horizonte infinito  % 
[Klqr,Pc]=lqr(A,B,Q,R); % Solução do problema LQR de horizonte Infinito
syslqr=ss(A-B*Klqr,B*0,C,D); % Sistema em malha fechada com LQR
[x2,t2] = initial(syslqr,x0,T); % Resposta do sistema de condição inicial
%%
figure, 
subplot(211), hold on, plot(t',Cx*x'),  plot(t2',Cx*x2','--','linewidth',2), xlabel('t (s)'), ylabel('x_i (cm)'), legend('x1(fin)','x2(fin)','x3(fin)','x1(inf)','x2(inf)','x3(inf)'), grid on
subplot(212), hold on, plot(t',Cv*x'),  plot(t2',Cv*x2','--','linewidth',2), xlabel('t (s)'), ylabel('v_i (cm/s)'), legend('v1(fin)','v2(fin)','v3(fin)','v1(inf)','v2(inf)','v3(inf)'), grid on
%% Reconstrução da entrada e do custo ótimo (aprox. via Euler direto)
u=zeros(2,length(t));J=0;
for i=1:length(t)
    for k=1:6*6
        vPk(k) = interp1(tvP(end:-1:1)',vP(end:-1:1,k)',t(i));
    end
    P=reshape(vPk,6,6);
    u(:,i)=-R\B'*P*x(i,:)';
    if i>1
        J=[J;J(i-1)+(t(i)-t(i-1))*(x(i,:)*Q*x(i,:).'+u(:,i)'*R*u(:,i))];
    end
end
J = J+x(end,:)*PT*x(end,:).';
u2=-Klqr*x2';
J2=0;
for i=2:length(t2)
    J2=[J2;J2(i-1)+(t2(i)-t2(i-1))*(x2(i,:)*Q*x2(i,:).'+u2(:,i)'*R*u2(:,i))];
end
figure, hold on, plot(t,u),  plot(t2,u2,'--','linewidth',2), xlabel('t (s)'), ylabel('u (N)'), legend('u1(fin)','u2(fin)','u1(inf)','u2(inf)'), grid on
figure, hold on, plot(t,J),  plot(t2,J2,'--','linewidth',2), xlabel('t (s)'), ylabel('J'), legend('J(fin)','J(inf)'), grid on