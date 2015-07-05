clear all;close all;clc;
%% Set model parameters
% A simple movement in 3-d space
m = [0, 0, 0, 1, 2, 3]';%[x, y, z, v_x, v_y, v_z]
n = length(m);
P = eye(n);
deltat = 1;
F = [eye(n/2),eye(n/2)*deltat;zeros(n/2),eye(n/2)];
Q = eye(n);
R = eye(n);
%% Set filter parameters
alpha = 1;
beta = 1;
kappa = 1;
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
ySE = 0;
mSE = 0;
for t = 0:100
%% Fake input
    y = [t, 2*t, 3*t, 1, 2, 3]'+chol(R)*randn(n,1);
    plot3(y(1),y(2),y(3),'rx');
    hold on;
    ySE = ySE + sum((y-[t, 2*t, 3*t, 1, 2, 3]').^2);
%% Prediction
    X = repmat(m,1,2*n+1)+sqrt(c)*[zeros(size(m)),chol(P),-chol(P)];
    Xhat = F*X;
    mminus = Xhat*W_m;
    Pminus = Xhat*W*Xhat' + Q;
%% Update
    Xminus = repmat(mminus,1,2*n+1)+sqrt(c)*[zeros(size(m)),chol(Pminus),-chol(Pminus)];
    Yminus = Xminus;
    mu = Yminus*W_m;
    S = Yminus*W*Yminus' + R;
    C = Xminus*W*Yminus';
    K = C/S;
    m = mminus + K*(y-mu);
    P = Pminus - K*S*K';
    plot3(m(1),m(2),m(3),'b.');
    hold on;
    mSE = mSE + sum((m-[t, 2*t, 3*t, 1, 2, 3]').^2);
end






