%% RBE 502 Prog Assignment 3
% Upasana Mahanti

clear; clc; close all;

t0 = 0;
tf = 10;
syms a0 a1 a2 a3 b0 b1 b2 b3 t x1 x2 x3 x4 A B lambda K

%% Part (e) :  Generation of desired trajectories 
D = [1 t0 t0^2 t0^3;0 1 2*t0 3*t0^2;1 tf tf^2 tf^3;0 1 2*tf 3*tf^2];
A1 = [a0;a1;a2;a3];  
A2 = [deg2rad(180);0;deg2rad(0);0];

B1 = [b0;b1;b2;b3];
B2 = [deg2rad(90);0;deg2rad(0);0];

A1 = inv(D)*A2;
B1 = inv(D)*B2;
   
qd=[A1';B1']*[1;t;t^2;t^3];
qd_dot=[A1';B1']*[0;1;2*t;3*t^2];
qd_ddot=[A1';B1']*[0;0;2;6*t];

[t,y] = ode45(@ode,[0 10],[deg2rad(200),deg2rad(125),0,0]); % angles in radian

for i=1:length(t)
    [~,u(i,:)] = ode(t(i,1),y(i,:));
end

% Plots
% State Trajectories

figure();
subplot(2,2,1)
hold on;
fplot(qd(1), [0 10])
plot(t,y(:,1));
xlabel('time t in sec');
ylabel('theta1');
hold off;

subplot(2,2,2)
hold on;
fplot(qd(2), [0 10])
plot(t,y(:,2));
xlabel('time t in sec');
ylabel('theta2');
hold off;

subplot(2,2,3)
hold on;
fplot(qd_dot(1), [0 10])
plot(t,y(:,3));
xlabel('time t in sec');
ylabel('theta1dot');
hold off;

subplot(2,2,4)
hold on;
fplot(qd_dot(2), [0 10])
plot(t,y(:,4));
xlabel('time t in sec');
ylabel('theta2dot');
hold off;



% Control input trajectories
figure()

subplot(2,1,1)
plot(t,u(:,1))
xlabel('time t in sec');
ylabel('u1');

subplot(2,1,2)
plot(t,u(:,2));
xlabel('time t in sec');
ylabel('u2');

