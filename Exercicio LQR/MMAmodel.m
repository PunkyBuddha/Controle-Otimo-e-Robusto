function dxdt = MMAmodel(t,x,A,B,R,tvP,vP)
    % Closed-loop model when controlled via finite-time LQR
    % P(t) is sent as a vector
    n=length(A);
    vPk=zeros(n*n,1); 
    for k=1:n*n
        % Interpolates P for the desired time t
        vPk(k) = interp1(tvP(end:-1:1)',vP(end:-1:1,k)',t);
    end
    P=reshape(vPk,n,n); %reconstructing the P(t) matrix
    dxdt = (A-B/R*B'*P)*x; %Closed loop system
end