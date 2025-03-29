clear; close all; clc;
try
    fclose(instrfindall);
catch
end
rosshutdown;

%% Timers
T_exp = 30; % Tempo de experimento
t_exp = tic;
T_run = 1/30; % Período do experimento
t_run = tic;
T_draw=0;
tempo = [];

%% Vetores de armazenamento
P = [0;0;0];
V = [0;0;0];
V_dot = [0;0;0];
nu = [0;0];
psi = 0;
pd = []; % Posição desejada para plot
pr = []; % Posição realizada para plot
pveld = []; % Velocidade desejada para plot
pvelr = []; % Velocidade realizada para plot
pu = []; % Esforço de controlador
er = []; % Erro para plot
ppsid = []; % Orientação desejada
ppsir = []; % Orientação realizada

%% Ganhos / Parametros
w = (2*pi)/15;
Kd = diag([3 6]);
Kp = diag([4.5 8.5]);
Ku = diag([.88 .88]);
Kv = diag([0.18227 0.17095]);
Kz = 1;
K_psi = 1;

%% Escolha do controlador
% 1 - LQR
flag_u = 2;

%% Modelo em Espaço de Estados
A = [0 0     1       0;
     0 0     0       1;
     0 0 -Kv(1,1)    0;
     0 0     0   -Kv(2,2)];
B = [0 0 Ku(1,1)   0;
     0 0    0   Ku(2,2)]';
C = eye(4);
D = zeros(4,2);

sysc = ss(A,B,C,D); % Sistema em tempo continuo
sysd = c2d(sysc,T_run,'tustin'); % Sistema em tempo discreto

%% Matrizes de ponderação LQR
% Regra de Bryson
Qii = [.01 .01 10 10].^-2;
Rii = [5 5].^2;

Q = diag(Qii);
R = diag(Rii);

Klqr = lqr(A,B,Q,R); % Solução do problema LQR de horizonte Infinito

%% Filtro de Kalman
% Inicialização
X = [P(1:2); V(1:2)];
P0 = diag([1 1 1 1].^2);
xkkm1=X;
Pkkm1=P0;

%% Matrizes de ponderação KF
% Regra de Bryson
Rw = diag([100 100 100 100].^2);
Rv = diag([100 100 100 100]).^2;

figure();

% while toc(t_exp) < T_exp
%     if toc(t_run) > T_run
%         tempo = [tempo toc(t_exp)];
%         dt = toc(t_run);
%         t_run = tic;
%         t = toc(t_exp);
%         t_corpo = tic;

for t = 0:T_run:T_exp
    tempo = [tempo t];
    t_run = tic;

    % Inovação (Kalman)
    X = [P(1:2); V(1:2)]; % Estados medidos em condições ideais
    y = C*X + D*nu; % Saída medida
    etak = y - C*xkkm1; 
    Sk = C*Pkkm1*C' + Rv;

    % Correção
    K = Pkkm1*C'/Sk;
    xkk = xkkm1 + K*etak;
    Pkk = (eye(4) - K*C)*Pkkm1;

    % Predição
    xkp1k = A*xkk;
    Pkp1k = A*Pkk*A' + Rw;

    % Volta ao Passo 1
    xkkm1 = xkp1k;
    Pkkm1 = Pkp1k;

    %% PLANEJADOR DE MOVIMENTO

    % % Lemniscata
    Pd = [sin(w*t); sin(2*w*t); 1]; % Posição desejada
    Vd = [cos(w*t)*w; cos(2*w*t)*2*w; 0]; % Velocidade desejada
    Vd_dot = [-sin(w*t)*w^2; -sin(2*w*t)*4*w^2; 0]; % Aceleração desejada

    Xr = [Pd(1:2); Vd(1:2)]; % Vetor de Estados de referência

    % % Orientação
    psid = [atan2(Vd(2),Vd(1)); 0]; % Orientação desejada

    pd = [pd Pd(1:3)]; % Armazenamento da posição desejada
    pveld = [pveld Vd(1:3)]; % Armazenamento da velocidade desejada
    ppsid = [ppsid psid(1)]; % Armazenamento da orientação desejada

    %% LEI DE CONTROLE

    if flag_u == 1
        nu = - Klqr*(X-Xr);
    else if flag_u == 2
            X = xkk;
            nu = - Klqr*(X-Xr);
    end
    end

    %% Controle em z

    Z_dot_ref = Vd(3) + Kz*(Pd(3) - P(3)); % Velocidade de referência em z

    Z_dot_ref = min(max(Z_dot_ref,-1),1); % Saturação de +-1m/s em z
    V(3) =  Z_dot_ref;

    %% Controle de orientação

    psi_til = psid(1) - psi; % Erro em psi

    % Filtro para otimização de caminho para orientação
    if abs(psi_til) > pi
        psi_til = psi_til - sign(psi_til)*2*pi;
    end

    psi_dot_ref = psid(2) + K_psi*psi_til; % Velocidade de referência em psi

    psi_dot_ref = min(max(psi_dot_ref,-1),1); % Limitador do controlador em psi

    u = [nu(1); nu(2); Z_dot_ref; psi_dot_ref]; % Vetor de comandos de controle Linear

    %% Simulação por discretização (x e y)
    Xkss = sysd.A*X + sysd.B*nu;
    P(1) = Xkss(1); P(2) = Xkss(2); V(1) = Xkss(3); V(2) = Xkss(4);

    %% Simulação por integração numérica (z e psi)
    P(3) = P(3) + V(3)*T_run; % Integracao de vel para achar pos em z
    psi = psi + psi_dot_ref*T_run; % Integração de psi_dot para encontrar orientação

    pr = [pr P(1:3)]; % Recebe a posição para plot de trajetoria
    pvelr = [pvelr V(1:3)];
    ppsir = [ppsir psi]; % Recebe orientação para calculo de erro

    pu = [pu u];

    %% Plot de trajetória online
    l = 0.1;
    T_draw=T_draw+T_run;
    if T_draw>0.5
        plot3(pd(1,:),pd(2,:),pd(3,:),'b--','LineWidth',1);
        hold on
        grid on
        axis([-2 2 -2 2 0 2])
        plot3(pr(1,:),pr(2,:),pr(3,:),'g','LineWidth',1);
        plot3(P(1),P(2),P(3),'g*','LineWidth',1);
        plot3([P(1), P(1) + l*cos(psi)], [X(2), X(2) + l*sin(psi)], [P(3) P(3)],'-r','LineWidth',1);
        plot3(Pd(1),Pd(2),Pd(3),'b*','LineWidth',1);
        T_draw=0;
        hold off
        drawnow
    end
end
% end

er = pd - pr; % Calculo de erros de posicionamento
err = ppsid - ppsir;
d_ec=sqrt(er(1,:).^2+er(2,:).^2); % Calculo de erro euclidiano

figure('Name','Graficos de erro')
subplot(4,1,1)
plot (tempo,er(1,:),'b','LineWidth',1);
xlabel('Tempo(s)');ylabel('Erro em x(m)');
grid on
subplot(4,1,2)
plot (tempo,er(2,:),'r','LineWidth',1)
xlabel('Tempo(s)');ylabel('Erro em y(m)');
grid on
subplot(4,1,3)
plot (tempo,er(3,:),'k','LineWidth',1)
xlabel('Tempo(s)');ylabel('Erro em z(m)');
grid on
subplot(4,1,4)
plot (tempo,err,'c','LineWidth',1)
xlabel('Tempo(s)');ylabel('Erro em psi(rad)');
grid on

figure('Name','Gráficos de posição vs Tempo')
subplot(4,1,1)
plot (tempo,pr(1,:),'r','LineWidth',1);
hold on
plot (tempo,pd(1,:),'b--','LineWidth',1);
xlabel('Tempo(s)');ylabel('Posição em x (m)');
legend('Posição realizada','Posição desejada');
grid on
subplot(4,1,2)
plot (tempo,pr(2,:),'r','LineWidth',1);
hold on
plot (tempo,pd(2,:),'b--','LineWidth',1);
xlabel('Tempo(s)');ylabel('Posição em y (m)');
legend('Posição realizada','Posição desejada');
grid on
subplot(4,1,3)
plot (tempo,pr(3,:),'r','LineWidth',1);
hold on
plot (tempo,pd(3,:),'b--','LineWidth',1);
xlabel('Tempo(s)');ylabel('Posição em z (m)');
legend('Posição realizada','Posição desejada');
grid on
subplot(4,1,4)
plot (tempo,ppsir,'r','LineWidth',1);
hold on
plot (tempo,ppsid,'b--','LineWidth',1);
xlabel('Tempo(s)');ylabel('Posição em psi (rad)');
legend('Posição realizada','Posição desejada');
grid on

figure('Name','Gráficos de velocidade vs Tempo')
subplot(3,1,1)
plot (tempo,pvelr(1,:),'r','LineWidth',1);
hold on
plot (tempo,pveld(1,:),'b--','LineWidth',1);
xlabel('Tempo(s)');ylabel('Velocidade em x (m/s)');
legend('Velocidade realizada','Velocidade desejada');
grid on
subplot(3,1,2)
plot (tempo,pvelr(2,:),'r','LineWidth',1);
hold on
plot (tempo,pveld(2,:),'b--','LineWidth',1);
xlabel('Tempo(s)');ylabel('Velocidade em y (m/s)');
legend('Velocidade realizada','Velocidade desejada');
grid on
subplot(3,1,3)
plot (tempo,pvelr(3,:),'r','LineWidth',1);
hold on
plot (tempo,pveld(3,:),'b--','LineWidth',1);
xlabel('Tempo(s)');ylabel('Velocidade em z (m/s)');
legend('Velocidade realizada','Velocidade desejada');
grid on

figure('Name','Esforço do controlador')
subplot(4,1,1)
plot (tempo,pu(1,:),'b');
xlabel('Tempo(s)');ylabel('Esforço do controlador em theta');
grid on
subplot(4,1,2)
plot (tempo,pu(2,:),'b');
xlabel('Tempo(s)');ylabel('Esforço do controlador em phi');
grid on
subplot(4,1,3)
plot (tempo,pu(3,:),'b');
xlabel('Tempo(s)');ylabel('Esforço do controlador em z');
grid on
subplot(4,1,4)
plot (tempo,pu(4,:),'b');
xlabel('Tempo(s)');ylabel('Esforço do controlador em psi');
grid on

