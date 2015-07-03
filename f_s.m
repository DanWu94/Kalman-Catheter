function [ sapo_k ] = f_s( s_km1 )
%Forward function
%% Get the catheter and plane state vectors
scat_km1 = s_km1.scat;
spi_km1 = s_km1.spi;
%% Predict the state of the plane. Primes indicate predictions.
sapopi_k = F_k*spi_km1;
%%
p_km1 = scat_km1.p;
v_km1 = scat_km1.v;
%%
uhatapo_k = sapopi_k.u
end

