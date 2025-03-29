clear all; close all; clc;
%% Modelo cinemático do avião
% v_{x,k+1} = v_{x,k}+w_{vx,k}
% p_{x,k+1} = p_{x,k}+T_s*v_{x,k}+w_{px,k};
% v_{y,k+1} = v_{y,k}+w_{vy,k}
% p_{y,k+1} = p_{y,k}+T_s*v_{y,k}+w_{py,k};
Ts=4;
A=[1 0 Ts 0 ;...
   0 1 0  Ts;...
   0 0 1  0 ;...
   0 0 0  1 ];
C=[eye(2) zeros(2)];
Rw=diag([1e-6 1e-6 20 20].^2);
Rv=diag([800 800]).^2;
% load reta.mat
load manobra.mat
y=x_m;
t=(0:(length(y)-1))*Ts;
x0=[-30000 -30000 0 0]';
P0=diag([100000 100000 100000 100000].^2);
% Passo 0: inicialização
xkkm1=x0;
Pkkm1=P0;
for k=1:length(y)
%Passo 1: Inovação
    etak = y(:,k)-C*xkkm1;
    Sk = C*Pkkm1*C'+Rv;
%Passo 2: Correção
    K{k} = Pkkm1*C'/Sk;
    xkk{k} = xkkm1+K{k}*etak;
    Pkk{k} = (eye(4)-K{k}*C)*Pkkm1;
%Passo 3: Predição
    xkp1k{k} = A*xkk{k};
    Pkp1k{k} = A*Pkk{k}*A'+Rw;
% Volta ao Passo 1
    xkkm1=xkp1k{k};
    Pkkm1=Pkp1k{k};
end
p=C*cell2mat(xkk);
figure, plot(x_c(:,1),x_c(:,2),'xg');
hold on, plot(x_m(1,:),x_m(2,:),'r*')
hold on, plot(p(1,:),p(2,:),'bo')
P=dare(A',C',Rw,Rv);
L=A*P*C'/(C*P*C'+Rv);
syskf=ss((A-L*C),L,C,0,Ts);
[xh,t]=lsim(syskf,x_m,t,x0);
hold on, plot(xh(:,1),xh(:,2),'sk')
%%
sysc=d2c(ss(A,[],C,0,Ts));
Ac=sysc.A;Bc=sysc.B;Cc=sysc.C;Dc=sysc.D;
Pc=care(Ac',Cc',Rw,Rv);
Lc=Pc*Cc'/Rv;
syskfc=ss(Ac-Lc*Cc,Lc,Cc,[]);
xhc=lsim(syskfc,x_m,t,x0);
hold on, plot(xhc(:,1),xhc(:,2),'k--')
figure, plot(t,sum(sqrt((x_c-x_m').^2),2),'g')
hold on, plot(t,sum(sqrt((x_c-p').^2),2),'b--')
hold on, plot(t,sum(sqrt((x_c-xh).^2),2),'r:')