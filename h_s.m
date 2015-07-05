function z_k = h_s( s_k )
%Measurement function
%%
O_k = s_k(1:3);
alpha_k = s_k(7);
beta_k = s_k(8);
gama_k = s_k(9);
p_k = s_k(13:15);
v_k = s_k(16:18);
%% Compute 1st plane vector
uhat_k = R_euler(alpha_k, beta_k, gama_k)*[1, 0, 0]';
%% Compute 2nd plane vector
vhat_k = R_euler(alpha_k, beta_k, gama_k)*[0, 1, 0]';
%% Compute intersection parameters
[lamda, u, v] = f_inter(p_k, v_k, uhat_k, vhat_k, O_k);
%%
z_k = zeros(14,1);
z_k(1:3) = O_k;
z_k(4) = alpha_k;
z_k(5) = beta_k;
z_k(6) = gama_k;
z_k(7:9) = p_k;
z_k(10:12) = v_k;
z_k(13) = u;
z_k(14) = v;
end

