function dPdt = riccatiODE(t, P1, A, C, Rv, Q)
P1 = reshape(P1, size(A)); %Convert from "n^2"-by-1 to "n"-by-"n"
dPdt = A*P1 + P1*A' - P1*C'*\Rv*C*P1 + Q; %Determine derivative
dPdt = dPdt(:); %Convert from "n"-by-"n" to "n^2"-by-1
end