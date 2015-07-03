function [lamda, u, v] = f_inter( p_k, v_k, uhat_k, vhat_k, O_k )
%intersection function
p_0 = p_k - O_k;
u = det([p_0,vhat_k,v_k])/det([uhat_k,vhat_k,v_k]);
v = det([uhat_k,p_0,v_k])/det([uhat_k,vhat_k,v_k]);
lamda = det([uhat_k,vhat_k,p_0])/det([uhat_k,vhat_k,v_k]);
end

