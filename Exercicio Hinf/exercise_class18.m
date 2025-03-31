clear all; close all; clc;
%% Modelo massa-mola-amortecedor triplo
m1=0.5;m2=0.5;m3=0.5;
k1=1;k2=2;
b1=0.02;b2=0.01;
AG=[  0         0        0      1         0        0   ; ...
      0         0        0      0         1        0   ; ...
      0         0        0      0         0        1   ; ...
   -k1/m1    +k1/m1      0   -b1/m1    +b1/m1      0   ; ...
   +k1/m2 -(k1+k2)/m2 +k2/m2 +b1/m2 -(b1+b2)/m2 +b2/m2 ; ...
      0      +k2/m3   -k2/m3    0      +b2/m3   -b2/m3   ...
];
BG=[   0     0   ;...
      0     0   ;...
      0     0   ;...
    +1/m1   0   ;...
      0     0   ;...
      0   -1/m3  ...
 ];
CG=[eye(3) zeros(3)];
DG=zeros(3,2);
%% Matrizes de ponderação do LQR
% Regra de Bryson
Qii=[2;2;2;1;1;1].^-2;
Rii=[20 20].^-2;
Q=diag(Qii);
R=diag(Rii);
%% Ruidos:
% Ruido de Processo (disturbio)
Rd = mdiag(0.01*eye(3), 0.001*eye(3));
% Ruido de medição
Rn = 0.1*eye(3);
%% Planta Generalizada
% Problema generalizado:
%   xdot = AG*x + BGu + \sqrt{R_d}*d
%   y = CG*x + DG*u + \sqrt{R_n}*n
%   z1 = \sqrt{Q}*x
%   z2 = \sqrt{R}*u
% z = [z1;z2]; w = [d;n];
sqrtQ = chol(Q);sqrtR = chol(R);sqrtRd = chol(Rd);sqrtRn = chol(Rn);
nx=size(AG,1); nu=size(BG,2); ny=size(CG,1);
nz1 = nx; nz2 = nu; nd = nx; nn=ny; nz = nz1+nz2; nw = nd+nn;
A = AG;
B1 = [sqrtRd zeros(nx,nn)];
B2=BG;
C1=[sqrtQ;zeros(nz2,nx)];
D11 = zeros(nz,nw);
D12 = [zeros(nz1,nu);sqrtR];
C2=CG;
D21 = [zeros(ny,nd) sqrtRn];
D22 = DG; % = zeros(ny,nu)

%% Verificação dos requisitos H2
Cc = ctrb(A,B2); % rank(Cc) = 6;
Oo = obsv(A,C2); % rank(Oo) = 6;
R1 = D12'*D12; % eig(R1) = [0.025;0.025]: R1=R1', lli(R1)>0 => R1>0
R2 = D21*D21'; % eig(R2) = [0.1;0.1;0.1]: R2=R2', lli(R2)>0 => R2>0

%% Solução LQG 
[Klqr,Pc]=lqr(AG,BG,Q,R); % Solução do problema LQR de horizonte Infinito
[LKF,Po]=lqe(AG,eye(6),CG,Rd,Rn); % Solução do problema KF estacionário
%% Solução H2
% Equação de Riccati 1:
%  AR1'*S1 + S1*AR1 - S1*BR1*S1 + CR1 = 0
AR1=A-B2/R1*D12'*C1;
BR1=B2/R1*B2';
CR1=C1'*C1-C1'*D12/R1*D12'*C1;
S1=are(AR1,BR1,CR1);
% Ganho de controle H2
FH2=R1\(B2'*S1 + D12'*C1);
% Equação de Riccati 2:
%  (AR2')'*S2 + S2*(AR2') - S2*BR2*S2 + CR2 = 0
AR2=A-B1*D21'/R2*C2;
BR2=C2'/R2*C2;
CR2=B1*B1'-B1*D21'/R2*D21*B1';
S2=are(AR2',BR2,CR2);
% Ganho de controle H2
LH2=(S2*C2' + B1*D21')/R2;
%% Solução H2 via matlab
% Planta generalizada:
AP = A;
BP = [B1 B2];
CP = [C1;C2];
DP = [D11 D12; D21 D22];
Ps = ss(AP,BP,CP,DP);
[Ks, Twz, nH2otima] = h2syn(Ps,ny,nu);
FH2matlab = -Ks.C;
LH2matlab = Ks.B;
%% Checagem norma H2 ótima
kp = [0.9:0.01:0.99 0.99:0.0001:1.01 1.01:0.01:1.1]; %Ganho proporcional ao ganho de controle H2
nH2 = 0*kp; % inicialização do vetor
for k=1:length(kp)
    nH2(k) = norm(lft(Ps,kp(k)*Ks),2); %calculo da norma H2 de Twz
end
figure, plot(kp,nH2), hold on, plot([kp(1) kp(end)], nH2otima*[1 1], '--k')
        xlabel('Ganho proporcional'), ylabel('||T_{wz}||_2'), xlim([0.9 1.1])
%% Solução Hinf 
%% Verificação dos requisitos
Cc = ctrb(A,B2); % rank(Cc) = 6;
Oo = obsv(A,C2); % rank(Oo) = 6;
ver1 = C1'*D12;  % C1'*D12 = zeros(nx,nu)
ver2 = B1*D21';  %  = zeros(nx,ny)
ver3 = D12'*D12; %  = 0.025*eye(2) => D12*u = [D12/norm(D12,2)]*[norm(D12,2)*u] = D12l*u
ver4 = D21*D21'; %  = 0.1*eye(3) => D21*w = [D21/norm(D21,2)]*[norm(D21,2)*w] = D21l*wl
%% Normaliozação da plata para D12'*D12=I e D21*D21'=I
%        xdot =          A*x +  (B1/n2D21)*(n2D21*w) +          B2*u
%   (z/n2D12) = (C1/n2D12)*x +           0*(n2D21*w) + (D12/n2D12)*u
%           y =         C2*x + (D21/n2D21)*(n2D21*w) +           0*u
% Controle: u = Ksinf*y; 
n2D21=norm(D21,2);n2D12=norm(D12,2);
An=A;         B1n=B1/n2D21;   B2n=B2;
C1n=C1/n2D12; D11n=D11;       D12n=D12/n2D12;
C2n=C2;       D21n=D21/n2D21; D22n=D22;
%% Método da Bisseção para norma ótima Hinf da malha fechada normalizada
gamma0=1e20; tol=0.001; kmax = 2000;
alpha = gamma0; beta = 1e-10; k=0;
while k<kmax && (alpha-beta)/beta > tol
    gamma = (alpha+beta)/2;
    % Equação de Riccati 1:
    %  AR1'*S1 + S1*AR1 - S1*BR1*S1 + CR1 = 0
    AR1i=An;
    BR1i=B2n*B2n' - gamma^(-2)*(B1n*B1n');
    CR1i=C1n'*C1n;
    S1i=are_mod(AR1i,BR1i,CR1i);
    % Equação de Riccati 2:
    %  (AR2')'*S2 + S2*(AR2') - S2*BR2*S2 + CR2 = 0
    AR2i=An;
    BR2i=C2n'*C2n - gamma^(-2)*(C1n'*C1n);
    CR2i=B1n*B1n';
    S2i=are_mod(AR2i',BR2i,CR2i); 
    check1 = ~isempty(S1i);
    check2 = ~isempty(S2i);
    if check1 && check2
        check3 = (max(real(eig(AR1i-BR1i*S1i))) < -1e-6); 
        check4 = (max(real(eig(AR2i-S2i*BR2i))) < -1e-6);
        check5 = (max(eig(S1i*S2i)) < gamma^2);
    else
        check3=0;check4=0;check5=0;
    end
    if check1 && check2 && check3 && check4 && check5
        alpha = gamma;
    else
        beta = gamma;
    end
    k = k+1;
end
gamma=alpha; % Norma Hinf ótima de T_{wy}
%% Controle ótimo Hinf
% Equação de Riccati 1:
%  AR1*S1 + S1*AR1' - S1*BR1*S1 + CR1 = 0
AR1i=An;
BR1i=B2n*B2n' - gamma^(-2)*(B1n*B1n');
CR1i=C1n'*C1n;
S1i=are_mod(AR1i,BR1i,CR1i);
% Equação de Riccati 2:
%  (AR2')*S2 + S2*(AR2')' - S2*BR2*S2 + CR2 = 0
AR2i=An;
BR2i=C2n'*C2n - gamma^(-2)*(C1n'*C1n);
CR2i=B1n*B1n';
S2i=are_mod(AR2i',BR2i,CR2i); 
Finf=B2n'*S1i;
Linf=S2i*C2n';
Zinf= eye(6)/(eye(6)-gamma^(-2)*S2i*S1i);
Ainf = An + gamma^(-2)*(B1n*B1n')*S1i-B2n*Finf-Zinf*Linf*C2n;
Ksinf=minreal(ss(Ainf,Zinf*Linf,-Finf,0));
%% Norma ótimo Kinf para planta real
nHinfotima_bis = n2D12*gamma*n2D21; % ||T_{wz}|| = ||(T_{\bar{w}z}*n2D21)|| = ||T_{\bar{w}z}||*n2D21
%% Solução Hinf via matlab
[Ksinfmatlab, Twzinf, nHinfotima_mat] =hinfsyn(Ps,ny,nu);
%% Checagem norma Hinf ótima
kp = [0.7:0.01:0.99 0.99:0.0001:1.02 1.01:0.01:1.3]; %Ganho proporcional ao ganho de controle H2
nHinf = 0*kp; % inicialização do vetor
for k=1:length(kp)
    nHinf(k) = norm(lft(Ps,kp(k)*Ksinf),inf); %calculo da norma Hinf de Twz
end
figure, plot(kp,nHinf), hold on, plot([kp(1) kp(end)], nHinfotima_bis*[1 1], '--k')
        xlabel('Ganho proporcional'), ylabel('||T_{wz}||_\infty'), xlim([0.7 1.3])
return
%% Simulação
% Condição inicial
x0=[1 3 7 -1 0 1]';
% Simulink
open('simPgenH2Hinf.slx');
sim('simPgenH2Hinf.slx');