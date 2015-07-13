clear all;close all;clc;
load data.mat;
%% Set model parameters
% A simple movement in 3-d space
m = [0, 0, 0, 0, 0, 12, ...
    0, 0, 0, 0, 0, 0, ...
    0, 0, 0, 1, 1, 1]';
% [O_1, O_2, O_3, Odot_1, Odot_2, Odot_3, 
% euler_1, euler_2, euler_3, eulerdot_1, eulerdot_2, eulerdot_3,
% p_1, p_2, p_3, v_1, v_2, v_3]
deltat = 1/30;
n = length(m);
P = 1*eye(n);
phi_s = 0.1;%don't know
C = [eye(3)*deltat^3/3, eye(3)*deltat^2/2; eye(3)*deltat^2/2, eye(3)*deltat].*phi_s;
Qpi = [C, zeros(size(C)); zeros(size(C)), C];
sigma_p_square = 0.1;%from 0.1 to 4.0
Qcat = eye(6)*sigma_p_square;
Q = zeros(n);
Q(1:12,1:12) = Qpi;
Q(13:18,13:18) = Qcat;
sigma_m_square = 0.05;%0.05, 2.0, 4.0
R = eye(14)*sigma_m_square;
%% Set filter parameters
alpha = 10;%determines spread of sigma points, usually small(e.g. 1e-3)
beta = 2;%incorporate prior knowledge of distribution, 2 is optimal for Gauss
kappa = 0;%secondary scaling parameter, usually set to 0
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
    %Xminus = repmat(mminus,1,2*n+1)+sqrt(c)*[zeros(size(m)),chol(Pminus),-chol(Pminus)];
    Yminus = zeros(14,2*n+1);
    for i = 1:2*n+1
        Yminus(:,i) = h_s(Xhat(:,i));
    end
    mu = Yminus*W_m;
    S = Yminus*W*Yminus' + R;
    C = Xhat*W*Yminus';
    %C = Xminus*W*Yminus';
    K = C/S;
    m = mminus + K*(y-mu);
    P = Pminus - K*S*K';
    %[V,D] = eig(P);
    %eps = 0.01;
    %d = diag(D);
    %d(d<=0) = eps;
    %P = V*diag(d)*V';
    plot3(m(13),m(14),m(15),'go');
    hold on;
%% Calculate squared error
    mSE = mSE + sum((m(13:15)-gr(7:9)).^2);
    ySE = ySE + sum((y(7:9)-gr(7:9)).^2);
end
title(['mSE = ',num2str(mSE),' ySE = ',num2str(ySE)]);






