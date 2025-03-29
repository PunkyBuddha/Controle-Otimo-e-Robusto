function dvPdt = riccatiDELQG(t,vP,A,B,Q,R)
    % Formulates the Riccati's Differential Equation (RDE)
    % P must be reshaped as vector to be integrated using ode45
    n=length(A);
    P1 = reshape(vP,n,n);
    dPdt = -(A*P1*C'/(C*P1*C' + Rv);
    dvPdt=reshape(dPdt,n*n,1);
end