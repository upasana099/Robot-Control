%% Trajectory Generation
clc
clear all
close all
syms a01 a11 a21 a31 a02 a12 a22 a32 t 'real'
syms theta1_d theta2_d theta1_dot_d theta2_dot_d q q1 q2 dq1 dq2 
syms a0 a1 a2 a3 b0 b1 b2 b3 t x1 x2 x3 x4


% Cubic Polynomial Trajectory
% Desired Trajectories

t0 = 0;
tf = 10;

D = [1 t0 t0^2 t0^3;0 1 2*t0 3*t0^2;1 tf tf^2 tf^3;0 1 2*tf 3*tf^2];                                                                                                                                
A1 = [a0;a1;a2;a3];  
A2 = [deg2rad(180);0;deg2rad(0);0];

B1 = [b0;b1;b2;b3];
B2 = [deg2rad(90);0;deg2rad(0);0];

A1 = inv(D)*A2;
B1 = inv(D)*B2;
 
qd=[A1';B1']*[1;t;t^2;t^3]
qd_dot=[A1';B1']*[0;1;2*t;3*t^2]
qd_ddot=[A1';B1']*[0;0;2;6*t]

qd_1 = qd(1,:);
qd_2 = qd(2,:);

qd_dot_1 = qd_dot(1,:);
qd_dot_2 = qd_dot(2,:);

qd_ddot_1 = qd_ddot(1,:);
qd_ddot_2 = qd_ddot(2,:);


fprintf(" Desired Trajectory of q1 \n")
disp (qd_1)


fprintf(" Desired Trajectory of q2 \n")
disp (qd_2)


fprintf("Desired Trajectory of q1_dot \n")
disp(qd_dot_1)


fprintf(" Desired Trajectory of q2_dot \n")

disp(qd_dot_2)

fprintf(" Desired Trajectory of q1_ddot \n")
disp(qd_ddot_1)

fprintf(" Desired Trajectory of q2_ddot \n")
disp(qd_ddot_2)

