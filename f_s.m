function [ sapo_k ] = f_s( s_km1,deltat )
%Forward function
%% Get the catheter and plane state vectors
spi_km1 = s_km1(1:12);
scat_km1 = s_km1(13:18);
%% Predict the state of the plane. Primes indicate predictions.
sapopi_k = stateTrans(spi_km1, deltat);
%%
p_km1 = scat_km1(1:3);
v_km1 = scat_km1(4:6);
%%
R = R_euler(spi_km1(7), spi_km1(8), spi_km1(9));
uhatapo_k = R*[1, 0, 0]';
vhatapo_k = R*[0, 1, 0]';
Oapo_k = sapopi_k(1:3);
%% Compute intersection parameters
% [lamda, u, v] = f_inter( p_km1, v_km1, uhatapo_k, vhatapo_k, Oapo_k );
%% Predict intersection point
% papo_k = u*uhatapo_k + v*vhatapo_k + Oapo_k;
% vapo_k = papo_k - p_km1;
papo_k = p_km1 + v_km1;
vapo_k = v_km1;
%%
sapo_k = zeros(size(s_km1));
sapo_k(1:12)= sapopi_k;
sapo_k(13:15) = papo_k;
sapo_k(16:18) = vapo_k;
end

