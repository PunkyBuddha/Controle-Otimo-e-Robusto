function dvPdt = riccatiDE(t,vP,A,B,Q,R)
    % Formulates the Riccati's Differential Equation (RDE)
    % P must be reshaped as vector to be integrated using ode45
    n=length(A);
    P = reshape(vP,n,n);
    dPdt = -(A'*P+P*A-P*B/R*B'*P+Q);
    dvPdt=reshape(dPdt,n*n,1);
end