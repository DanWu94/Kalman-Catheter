clear all;close all;clc;
% Generate fake data for simulation
%% Set simulation parameters
O_0 = [0, 0, 0]';
Odot = [0, 0, 1]';
euler = [0, 0, 0]';
R = R_euler(euler(1), euler(2), euler(3));
uhat = R(:,1);
vhat = R(:,2);
sigma = 0.1;
deltat = 1;
T = 10;
%% Data generating
data = zeros(14,T/deltat);
for t = 1:T/deltat
    data(1:3,t) = O_0 + Odot.*t*deltat;
    data(4:6,t) = euler;
    data(13,t) = 0.01*t*t*t*deltat;
    data(14,t) = t*deltat;
    data(7:9,t) = data(1:3,t) + data(13,t)*uhat + data(14,t)*vhat;
    if t == 1
        data(10:12,t) = data(7:9,t);
    else
        data(10:12,t) = data(7:9,t)-data(7:9,t-1);
    end
end
obdata = data + sigma*randn(size(data));
figure;
for t = 1:T/deltat
   plot3(data(7,t),data(8,t),data(9,t),'b.');
   hold on;
   plot3(obdata(7,t),obdata(8,t),obdata(9,t),'rx');
   hold on;
end
save('data.mat','data','obdata');