function [ s_k ] = stateTrans( s_km1, deltat )
%State Transition
s_k.O = s_km1.O + s_km1.Ofd*deltat;
s_k.Ofd = s_km1.Ofd;
s_k.alpha = s_km1.alpha + s_km1.alphafd*deltat;
s_k.alphafd = s_km1.alphafd;
s_k.beta = s_km1.beta + s_km1.betafd*deltat;
s_k.betafd = s_km1.betafd;
s_k.gama = s_km1.gama + s_km1.gamafd*deltat;
s_k.gamafd = s_km1.gamafd;
end

