%% Programming Assignmenet 1
%% Upasana Mahanti
clc;
clear all;
close all;

%% 1.2 equations of motion for the RRBot

syms theta_1(t) theta_2(t) theta_ddot_1 theta_ddot_2 'real';
syms r1 r2 l1 l2 m1 m2 I1 I2 tau1 tau2 w1 w2 g KE PE 'real';

% variables 
x1 = r1*sin(theta_1);
y1 = r1*cos(theta_1);
x2 = l1*sin(theta_1) + r2*sin(theta_1 + theta_2);
y2 = l1*cos(theta_1) + r2*cos(theta_1 + theta_2);

xd1 = diff(x1,t);
yd1 = diff(y1,t);
xd2 = diff(x2,t);
yd2 = diff(y2,t);

q = sym('q',[2,1]);
q(1) = theta_1;
q(2) = theta_2;
disp(q);
q_dot = diff(q,t);
disp(q_dot);

u = sym('u',[2,1]);
u(1) = tau1;
u(2) = tau2;

v1 = sqrt(xd1^2 + yd1^2);
v2 = sqrt(xd2^2 + yd2^2);

w1 = jacobian(theta_1,t);
w2 = jacobian(theta_1+theta_2,t);
%disp(w1,w2);


% Kinetic Energy
KE = simplify((m1*(v1^2))/2 + (I1*w1^2)/2 + (m2*(v2^2))/2 + (I2*w2^2)/2);
disp(KE);

h1 = y1;
h2 = y2;

% Potential Energy
PE = m1*g*h1 + m2*g*h2;
disp(PE);

%Langrangian 
L = KE-PE;

dl_dtheta1 = jacobian(L,q(1));
dl_dtheta2 = jacobian(L,q(2));

dl_dtheta1_dot = jacobian(L, q_dot(1));
dl_dtheta2_dot = jacobian(L, q_dot(2));

ddl_dtheta1_dot_dt = jacobian(dl_dtheta1_dot,[q(1);q_dot(1)])*[q_dot(1);theta_ddot_1] + jacobian(dl_dtheta1_dot,[q(2);q_dot(2)])*[q_dot(2);theta_ddot_2];
ddl_dtheta2_dot_dt = jacobian(dl_dtheta2_dot,[q(1);q_dot(1)])*[q_dot(1);theta_ddot_1] + jacobian(dl_dtheta2_dot,[q(2);q_dot(2)])*[q_dot(2);theta_ddot_2];

eq1= ddl_dtheta1_dot_dt - dl_dtheta1 - tau1;

eq2 = ddl_dtheta2_dot_dt - dl_dtheta2 - tau2;

disp(eq1);
disp(eq2);

%% 1.2 (b) State-Space Representation

% Solving eq1 and eq2 to get theta_ddot_1 and theta_ddot_2

sol = solve([eq1==0, eq2==0], [theta_ddot_1, theta_ddot_2]);

display(sol.theta_ddot_1);
display(sol.theta_ddot_2);

%% 1.2 (c) ODE

[t, y] = ode45(@ode_planner_2link, [0, 10], [deg2rad(30), deg2rad(45), 0, 0]);


%% Trajectory Plots

figure('Name','Link 1 position');
%subplot(2,2,1);
plot(t,y(:,1));
xlabel('Time (in secs)');
ylabel('theta1 (in rads)');

figure('Name','Link 2 position');
%subplot(2,2,2);
plot(t,y(:,2));
xlabel('Time (in secs)');
ylabel('theta2 (in rads)');

figure('Name','Link 1 velocity');
%subplot(2,2,3);
plot(t,y(:,3));
xlabel('Time (in secs)');
ylabel('theta1 dot (in rads per sec)');

figure('Name','Link 2 velocity');
%subplot(2,2,4);
plot(t,y(:,4));
xlabel('Time (in secs)');
ylabel('theta2 dot (in rads per sec)')
