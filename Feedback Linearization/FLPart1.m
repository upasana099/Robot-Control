%% RBE 502 Prog Assignment 3
% Upasana Mahanti

clear; close; clc;

syms theta1 theta2 theta1_dot theta2_dot theta1_Ddot theta2_Ddot m1 m2 l1 l2 g I1 I2 r1 r2 u1 u2 'real'; 
syms k1 k2 Kp1 Kp2 Kd1 Kd2
syms theta1_d theta2_d theta1_dot_d theta2_dot_d
syms a0 a1 a2 a3 b0 b1 b2 b3 t x1 x2 x3 x4 A B lambda K


%% Part(a) : Generation of a cubic polynomial trajectory for the first and second joint of the robot
t0 = 0;
tf = 10;

D = [1 t0 t0^2 t0^3;0 1 2*t0 3*t0^2;1 tf tf^2 tf^3;0 1 2*tf 3*tf^2];
A1 = [a0;a1;a2;a3];  
A2 = [deg2rad(180);0;deg2rad(0);0];

B1 = [b0;b1;b2;b3];
B2 = [deg2rad(90);0;deg2rad(0);0];

A1 = inv(D)*A2
B1 = inv(D)*B2

% disp([A1';B1'])   

% Mathematical formulation and derivation:

qd=[A1';B1']*[1;t;t^2;t^3]
qd_dot=[A1';B1']*[0;1;2*t;3*t^2]
qd_ddot=[A1';B1']*[0;0;2;6*t]

% Desired Trajectories

figure
subplot(2,3,1)
fplot(qd(1), [0 10])
xlabel('time t in sec');
ylabel('theta1d');
subplot(2,3,2)
fplot(qd_dot(1), [0 10])
xlabel('time t in sec');
ylabel('theta1ddot');
subplot(2,3,3)
fplot(qd_ddot(1), [0 10])
xlabel('time t in sec');
ylabel('theta1dddot');
subplot(2,3,4)
fplot(qd(2), [0 10])
xlabel('time t in sec');
ylabel('theta2d');
subplot(2,3,5)
fplot(qd_dot(2), [0 10])
xlabel('time t in sec');
ylabel('theta2ddot');
subplot(2,3,6)
fplot(qd_ddot(2), [0 10])
xlabel('time t in sec');
ylabel('theta2dddot');


%% Part (b) : Manipulator Form

u = [u1;u2];
q=[theta1; theta2];
dq=[theta1_dot; theta2_dot];
ddq=[theta1_Ddot; theta2_Ddot]; 

% Considering the equations of motion derived for the robot in Programming Assignment 1

u(1)=theta1_Ddot*(I1 + I2 + (m1*(2*r1^2*cos(theta1)^2 + 2*r1^2*sin(theta1)^2))/2 + (m2*(2*(r2*cos(theta1 + theta2) + l1*cos(theta1))^2 + 2*(r2*sin(theta1 + theta2) + l1*sin(theta1))^2))/2) + theta2_Ddot*(I2 + (m2*(2*r2*sin(theta1 + theta2)*(r2*sin(theta1 + theta2) + l1*sin(theta1)) + 2*r2*cos(theta1 + theta2)*(r2*cos(theta1 + theta2) + l1*cos(theta1))))/2) - g*m2*(r2*sin(theta1 + theta2) + l1*sin(theta1)) - (m2*theta2_dot*(2*r2*sin(theta1 + theta2)*(l1*theta1_dot*cos(theta1) + r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot)) - 2*r2*cos(theta1 + theta2)*(l1*theta1_dot*sin(theta1) + r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot)) + 2*r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot)*(r2*cos(theta1 + theta2) + l1*cos(theta1)) - 2*r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot)*(r2*sin(theta1 + theta2) + l1*sin(theta1))))/2 - g*m1*r1*sin(theta1);

u(2)=theta2_Ddot*(I2 + (m2*(2*r2^2*cos(theta1 + theta2)^2 + 2*r2^2*sin(theta1 + theta2)^2))/2) + (m2*(2*r2*sin(theta1 + theta2)*(l1*theta1_dot*cos(theta1) + r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot))*(theta1_dot + theta2_dot) - 2*r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot)*(l1*theta1_dot*sin(theta1) + r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot))))/2 + theta1_Ddot*(I2 + (m2*(2*r2*sin(theta1 + theta2)*(r2*sin(theta1 + theta2) + l1*sin(theta1)) + 2*r2*cos(theta1 + theta2)*(r2*cos(theta1 + theta2) + l1*cos(theta1))))/2) - (m2*theta2_dot*(2*r2*sin(theta1 + theta2)*(l1*theta1_dot*cos(theta1) + r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot)) - 2*r2*cos(theta1 + theta2)*(l1*theta1_dot*sin(theta1) + r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot))))/2 - g*m2*r2*sin(theta1 + theta2);


g_q(1,1) = subs(u(1),[theta1_dot,theta2_dot,theta1_Ddot,theta2_Ddot],[0,0,0,0]);
g_q(2,1) = subs(u(2),[theta1_dot,theta2_dot,theta1_Ddot,theta2_Ddot],[0,0,0,0]);
g_q=simplify(g_q)

C_q_qdot1 = u(1)-g_q(1,1);
C_q_qdot(1,1) = subs(C_q_qdot1,[theta1_Ddot,theta2_Ddot],[0,0]);
C_q_qdot2 = u(2)-g_q(2,1);
C_q_qdot(2,1) = subs(C_q_qdot2,[theta1_Ddot,theta2_Ddot],[0,0]);
C_q_qdot=simplify(C_q_qdot)

M_q1 = u(1)- g_q(1,1);
M_q2 = u(2)- g_q(2,1);
M_q(1,1) = subs(M_q1,[theta1_dot,theta2_dot,theta2_Ddot,theta1_Ddot],[0,0,0,1]);
M_q(2,1) = subs(M_q2,[theta1_dot,theta2_dot,theta2_Ddot,theta1_Ddot],[0,0,0,1]);
M_q(1,2) = subs(M_q1,[theta1_dot,theta2_dot,theta1_Ddot,theta2_Ddot],[0,0,0,1]);
M_q(2,2) = subs(M_q2,[theta1_dot,theta2_dot,theta1_Ddot,theta2_Ddot],[0,0,0,1]);
M_q=simplify(M_q)

% Manipulator Form:
u = M_q*ddq + C_q_qdot + g_q

%% Part (c) : Symbolic feedback linearization of the robot

X=sym('X',[4,1]);
X(1) = theta1;
X(2) = theta2;
X(3) = theta1_dot; 
X(4) = theta2_dot;

Xd=[qd;qd_dot]

% Eigenvalue placement method

A=[0 1;0 0];
B=[0;1];
K=[Kp1 0 Kp2 0;0 Kd1 0 Kd2]
lambda1=[-3 -2];
K1 = place(A,B,lambda1);
lambda2=[-1 -5];
K2 = place(A,B,lambda2);
K = [diag(K1) diag(K2)]


% State-feedback control for the virtual control input

v = -K*(X-Xd)+qd_ddot
u = M_q*v + C_q_qdot + g_q
