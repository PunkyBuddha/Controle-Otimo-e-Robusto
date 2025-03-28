clear; close all; clc;
try
    fclose(instrfindall);
catch
end
rosshutdown;

% timers
T_exp = 120; % Tempo de experimento
t_exp = tic;
T_run = 1/30; % Período do experimento
t_run = tic;
T_draw=0;
tempo = [];

% Vetores de armazenamento
X = [0;0;0];
X_dot = [0;0;0];
X_2dot = [0;0;0];
psi = 0;
pd = []; % Posição desejada para plot
pr = []; % Posição realizada para plot
pveld = []; % Velocidade desejada para plot
pvelr = []; % Velocidade realizada para plot
paccr = []; % Aceleração realizada para plot
paccd = []; % Aceleração desejada para plot
pu = []; % Esforço de controlador
er = []; % Erro para plot
ppsid = []; % Orientação desejada
ppsir = []; % Orientação realizada

%% Ganhos / Parametros
w = (2*pi)/30;
Kd = diag([3 6]);
Kp = diag([4.5 8.5]);
Ku = diag([.88 .88 1]);
Kv = diag([0.18227 0.17095 1]);
Kz = 1;
K_psi = 1;

%% Modelo em Espaço de Estados
A = [0 0 0     1       0       0;
     0 0 0     0       1       0;
     0 0 0     0       0       1;
     0 0 0 -Kv(1,1)    0       0;
     0 0 0     0   -Kv(2,2)    0;
     0 0 0     0       0   -Kv(3,3)];
B = [0 0 0 Ku(1,1)   0      0;
     0 0 0    0   Ku(2,2)   0;
     0 0 0    0      0   Ku(3,3)]';
C = zeros(1,6);
D = zeros(1,3);

%% Matrizes de ponderação LQR
% Regra de Bryson
Qii = [10 10 10 5 5 5].^-2;
Rii = [1 1 1].^2;

Q = diag(Qii);
R = diag(Rii);

x0 = [0 0 0 0 0 0]';

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

    %% PLANEJADOR DE MOVIMENTO

    % % Lemniscata
    Xd = [sin(w*t); sin(2*w*t); 1]; % Posição desejada
    Xd_dot = [cos(w*t)*w; cos(2*w*t)*2*w; 0]; % Velocidade desejada
    Xd_2dot = [-sin(w*t)*w^2; -sin(2*w*t)*4*w^2; 0]; % Aceleração desejada

    % % Orientação
    psid = [0; 0]; % Orientação desejada

    pd = [pd Xd(1:3)]; % Armazenamento da posição desejada
    pveld = [pveld Xd_dot(1:3)]; % Armazenamento da velocidade desejada
    paccd = [paccd Xd_2dot(1:3)]; % Armazenamento da aceleração desejada
    ppsid = [ppsid psid(1)]; % Armazenamento da orientação desejada

    %% LEI DE CONTROLE

    R = [cos(psi) -sin(psi);
        sin(psi) cos(psi)];

    X_til = Xd - X; % Erro de posicionamento
    X_dot_til = Xd_dot - X_dot; % Erro de Velocidade

    x_2dot_ref = Xd_2dot(1:2) + Kd*X_dot_til(1:2) + Kp*X_til(1:2); % Aceleração de referência

    nu = inv(R*Ku)*(x_2dot_ref + Kv*X_dot(1:2)); % Lei de controle linear

    theta = min(max(nu(1),-1),1); % Saturação de +-1 em theta
    phi = min(max(nu(2),-1),1); % Saturação de +-1 em phi

    nu = [theta; phi];

    %% Controle em z

    Z_dot_ref = Xd_dot(3) + Kz*X_til(3); % Velocidade de referência em z

    Z_dot_ref = min(max(Z_dot_ref,-1),1); % Saturação de +-1m/s em z
    X_dot(3) =  Z_dot_ref;

    %% Controle de orientação

    psi_til = psid(1) - psi; % Erro em psi

    % Filtro para otimização de caminho para orientação
    if abs(psi_til) > pi
        psi_til = psi_til - sign(psi_til)*2*pi;
    end

    psi_dot_ref = psid(2) + K_psi*psi_til; % Velocidade de referência em psi

    psi_dot_ref = min(max(psi_dot_ref,-1),1); % Limitador do controlador em psi

    u = [theta; -phi; Z_dot_ref; psi_dot_ref]; % Vetor de comandos de controle Linear

    X_2dot(1:2) = R*Ku*nu - Kv*X_dot(1:2);

    


    X_2dot(1) = min(max(X_2dot(1),-1),1);
    X_2dot(2) = min(max(X_2dot(2),-1),1);

    %% SIMULACAO

    X_dot(1:2) = X_dot(1:2) + X_2dot(1:2)*T_run; % Integracao de acc para achar vel
    X(1:2) = X(1:2) + X_dot(1:2)*T_run; % Integracao de vel para achar pos
    X(3) = X(3) + X_dot(3)*T_run; % Integracao de vel para achar pos em z
    psi = psi + psi_dot_ref*T_run; % Integração de psi_dot para encontrar orientação

    pr = [pr X(1:3)]; % Recebe a posição para plot de trajetoria
    pvelr = [pvelr X_dot(1:3)];
    ppsir = [ppsir psi]; % Recebe orientação para calculo de erro

    pu = [pu u];
    paccr = [paccr X_2dot];

    l = 0.1;
    T_draw=T_draw+T_run;
    if T_draw>0.5
        plot3(pd(1,:),pd(2,:),pd(3,:),'b--','LineWidth',1);
        hold on
        grid on
        axis([-2 2 -2 2 0 2])
        plot3(pr(1,:),pr(2,:),pr(3,:),'g','LineWidth',1);
        plot3(X(1),X(2),X(3),'g*','LineWidth',1);
        plot3([X(1), X(1) + l*cos(psi)], [X(2), X(2) + l*sin(psi)], [X(3) X(3)],'-r','LineWidth',1);
        plot3(Xd(1),Xd(2),Xd(3),'b*','LineWidth',1);
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

% figure('Name','Gráficos de posição vs Tempo')
% subplot(4,1,1)
% plot (tempo,pr(1,:),'r','LineWidth',1);
% hold on
% plot (tempo,pd(1,:),'b--','LineWidth',1);
% xlabel('Tempo(s)');ylabel('Posição em x (m)');
% legend('Posição realizada','Posição desejada');
% grid on
% subplot(4,1,2)
% plot (tempo,pr(2,:),'r','LineWidth',1);
% hold on
% plot (tempo,pd(2,:),'b--','LineWidth',1);
% xlabel('Tempo(s)');ylabel('Posição em y (m)');
% legend('Posição realizada','Posição desejada');
% grid on
% subplot(4,1,3)
% plot (tempo,pr(3,:),'r','LineWidth',1);
% hold on
% plot (tempo,pd(3,:),'b--','LineWidth',1);
% xlabel('Tempo(s)');ylabel('Posição em z (m)');
% legend('Posição realizada','Posição desejada');
% grid on
% subplot(4,1,4)
% plot (tempo,ppsir,'r','LineWidth',1);
% hold on
% plot (tempo,ppsid,'b--','LineWidth',1);
% xlabel('Tempo(s)');ylabel('Posição em psi (rad)');
% legend('Posição realizada','Posição desejada');
% grid on
% 
% figure('Name','Gráficos de velocidade vs Tempo')
% subplot(3,1,1)
% plot (tempo,pvelr(1,:),'r','LineWidth',1);
% hold on
% plot (tempo,pveld(1,:),'b--','LineWidth',1);
% xlabel('Tempo(s)');ylabel('Velocidade em x (m/s)');
% legend('Velocidade realizada','Velocidade desejada');
% grid on
% subplot(3,1,2)
% plot (tempo,pvelr(2,:),'r','LineWidth',1);
% hold on
% plot (tempo,pveld(2,:),'b--','LineWidth',1);
% xlabel('Tempo(s)');ylabel('Velocidade em y (m/s)');
% legend('Velocidade realizada','Velocidade desejada');
% grid on
% subplot(3,1,3)
% plot (tempo,pvelr(3,:),'r','LineWidth',1);
% hold on
% plot (tempo,pveld(3,:),'b--','LineWidth',1);
% xlabel('Tempo(s)');ylabel('Velocidade em z (m/s)');
% legend('Velocidade realizada','Velocidade desejada');
% grid on
% 
% figure('Name','Gráficos de Aceleração vs Tempo')
% subplot(2,1,1)
% plot (tempo,paccr(1,:),'r','LineWidth',1);
% hold on
% plot (tempo,paccd(1,:),'b--','LineWidth',1);
% xlabel('Tempo(s)');ylabel('Aceleração em x (m/s^2)');
% legend('Aceleração realizada','Aceleração desejada');
% grid on
% subplot(2,1,2)
% plot (tempo,paccr(2,:),'r','LineWidth',1);
% hold on
% plot (tempo,paccd(2,:),'b--','LineWidth',1);
% xlabel('Tempo(s)');ylabel('Aceleração em y (m/s^2)');
% legend('Aceleração realizada','Aceleração desejada');
% grid on
% 
% figure('Name','Esforço do controlador')
% subplot(4,1,1)
% plot (tempo,pu(1,:),'b');
% xlabel('Tempo(s)');ylabel('Esforço do controlador em theta');
% grid on
% subplot(4,1,2)
% plot (tempo,pu(2,:),'b');
% xlabel('Tempo(s)');ylabel('Esforço do controlador em phi');
% grid on
% subplot(4,1,3)
% plot (tempo,pu(3,:),'b');
% xlabel('Tempo(s)');ylabel('Esforço do controlador em z');
% grid on
% subplot(4,1,4)
% plot (tempo,pu(4,:),'b');
% xlabel('Tempo(s)');ylabel('Esforço do controlador em psi');
% grid on

