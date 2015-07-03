function z_k = h_s( s_k )
%Measurement function
%%
p_k = s_k.scat.p;
v_k = s_k.scat.v;
O_k = s_k.spi.O;
alpha_k = s_k.spi.alpha;
beta_k = s_k.spi.beta;
gama_k = s_k.spi.gama;
%% Compute 1st plane vector
uhat_k = R_euler(alpha_k, beta_k, gama_k)*[1, 0, 0]';
%% Compute 2nd plane vector
vhat_k = R_euler(alpha_k, beta_k, gama_k)*[0, 1, 0]';
%% Compute intersection parameters
[lamda, u, v] = f_inter(p_k, v_k, uhat_k, vhat_k, O_k);
%%
z_k.u = u;
z_k.v = v;
z_k.p_k = p_k;
z_k.v_k = v_k;
z_k.O_k = O_k;
z_k.alpha_k = alpha_k;
z_k.beta_k = beta_k;
z_k.gama_k = gama_k;
end

