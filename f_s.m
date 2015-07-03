function [ sapo_k ] = f_s( s_km1,deltat )
%Forward function
%% Get the catheter and plane state vectors
scat_km1 = s_km1.scat;
spi_km1 = s_km1.spi;
%% Predict the state of the plane. Primes indicate predictions.
sapopi_k = stateTrans(spi_km1, deltat);
%%
p_km1 = scat_km1.p;
v_km1 = scat_km1.v;
%%
uhatapo_k = sapopi_k.u;
vhatapo_k = sapopi_k.v;
Oapo_k = sapopi_k.O;
%% Compute intersection parameters
[lamda, u, v] = f_inter( p_km1, v_km1, uhatapo_k, vhatapo_k, Oapo_k );
%% Predict intersection point
papo_k = u*uhatapo_k + v*vhatapo_k + Oapo_k;
vapo_k = papo_k - p_km1;
%%
sapocat_k.p = papo_k;
sapocat_k.v = vapo_k;
%%
sapo_k.scat = sapocat_k;
sapo_k.pi = sapopi_k;
end

