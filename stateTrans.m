function [ s_k ] = stateTrans( s_km1, deltat )
%State Transition
A = [eye(3),eye(3)*deltat;zeros(3),eye(3)];
F = [A, zeros(size(A));zeros(size(A)),A];
s_k = F*s_km1;
end

