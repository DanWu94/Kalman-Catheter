clear all;close all;clc;
load data.mat;
%% Set model parameters
% A simple movement in 3-d space
m = [0, 0, 0, 0, 0, 1, ...
    0, 0, 0, 0, 0, 0, ...
    0, 0, 0, 1, 1, 1]';
% [O_1, O_2, O_3, Odot_1, Odot_2, Odot_3, 
% euler_1, euler_2, euler_3, eulerdot_1, eulerdot_2, eulerdot_3,
% p_1, p_2, p_3, v_1, v_2, v_3]
deltat = 1;
n = length(m);
P = 1*eye(n);
phi_s = 1;
C = [eye(3)*deltat^3/3, eye(3)*deltat^2/2; eye(3)*deltat^2/2, eye(3)*deltat]*phi_s;
Qpi = [C, zeros(size(C)); zeros(size(C)), C];
sigma_p = 1;
Qcat = eye(6)*sigma_p^2;
Q = zeros(n);
Q(1:12,1:12) = Qpi;
Q(13:18,13:18) = Qcat;
sigma_m = 0.1;
R = eye(14)*sigma_m^2;
%% Set filter parameters
alpha = 1;
beta = 1;
kappa = 10;
lamda = alpha^2*(n+kappa)-n;
c = alpha^2*(n+kappa);
W_m = zeros(2*n+1,1);
W_c = zeros(2*n+1,1);
W_m(1) = lamda/(n+lamda);
W_c(1) = lamda/(n+lamda)+(1-alpha^2+beta);
for i = 1:2*n
   W_m(i+1) = 1/(2*(n+lamda));
   W_c(i+1) = 1/(2*(n+lamda));
end
W = (eye(2*n+1)-repmat(W_m,1,2*n+1))*diag(W_c)*(eye(2*n+1)-repmat(W_m,1,2*n+1))';
%% Interation
figure;
for t = 1:size(data,2)
%% Fake input
    y = obdata(:,t);    
    plot3(y(7),y(8),y(9),'rx');
    hold on;
%% Plot ground truth
    gr = data(:,t);
    plot3(gr(7),gr(8),gr(9),'b.');
    hold on;
%% Prediction
    X = repmat(m,1,2*n+1)+sqrt(c)*[zeros(size(m)),chol(P),-chol(P)];
    Xhat = zeros(size(X));
    for i = 1:2*n+1
        Xhat(:,i) = f_s(X(:,i),deltat);
    end
    mminus = Xhat*W_m;
    Pminus = Xhat*W*Xhat' + Q;
%% Update
    Xminus = repmat(mminus,1,2*n+1)+sqrt(c)*[zeros(size(m)),chol(Pminus),-chol(Pminus)];
    Yminus = zeros(14,2*n+1);
    for i = 1:2*n+1
        Yminus(:,i) = h_s(Xminus(:,i));
    end
    mu = Yminus*W_m;
    S = Yminus*W*Yminus' + R;
    C = Xminus*W*Yminus';
    K = C/S;
    m = mminus + K*(y-mu);
    P = Pminus - K*S*K';
    [V,D] = eig(P);
    eps = 0.01;
    d = diag(D);
    d(d<=0) = eps;
    P = V*diag(d)*V';
    plot3(m(13),m(14),m(15),'go');
    hold on;
end






