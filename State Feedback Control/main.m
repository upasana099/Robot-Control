%% Programming Assignmenet 1
%% Upasana Mahanti
clc;
clear all;
close all;

%%
syms theta1 theta2 theta1_dot theta2_dot theta1_ddot theta2_ddot m1 m2 l1 l2 g I1 I2 r1 r2 u1 u2 'real'; 
syms X1 X2 X3 X4 lambda

%% Q2.2 (a) Equilibrium Points

u = sym('u',[2,1]);
u(2) = u2;
u(1) = u1;


eq1 = theta1_ddot*(I1 + I2 + (m1*(2*r1^2*cos(theta1)^2 + 2*r1^2*sin(theta1)^2))/2 + (m2*(2*(r2*cos(theta1 + theta2) + l1*cos(theta1))^2 + 2*(r2*sin(theta1 + theta2) + l1*sin(theta1))^2))/2) - u1 + theta2_ddot*(I2 + (m2*(2*r2*sin(theta1 + theta2)*(r2*sin(theta1 + theta2) + l1*sin(theta1)) + 2*r2*cos(theta1 + theta2)*(r2*cos(theta1 + theta2) + l1*cos(theta1))))/2) - g*m2*(r2*sin(theta1 + theta2) + l1*sin(theta1)) - (m2*theta2_dot*(2*r2*sin(theta1 + theta2)*(l1*theta1_dot*cos(theta1) + r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot)) - 2*r2*cos(theta1 + theta2)*(l1*theta1_dot*sin(theta1) + r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot)) + 2*r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot)*(r2*cos(theta1 + theta2) + l1*cos(theta1)) - 2*r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot)*(r2*sin(theta1 + theta2) + l1*sin(theta1))))/2 - g*m1*r1*sin(theta1);
eq2 = theta2_ddot*(I2 + (m2*(2*r2^2*cos(theta1 + theta2)^2 + 2*r2^2*sin(theta1 + theta2)^2))/2) - u2 + (m2*(2*r2*sin(theta1 + theta2)*(l1*theta1_dot*cos(theta1) + r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot))*(theta1_dot + theta2_dot) - 2*r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot)*(l1*theta1_dot*sin(theta1) + r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot))))/2 + theta1_ddot*(I2 + (m2*(2*r2*sin(theta1 + theta2)*(r2*sin(theta1 + theta2) + l1*sin(theta1)) + 2*r2*cos(theta1 + theta2)*(r2*cos(theta1 + theta2) + l1*cos(theta1))))/2) - (m2*theta2_dot*(2*r2*sin(theta1 + theta2)*(l1*theta1_dot*cos(theta1) + r2*cos(theta1 + theta2)*(theta1_dot + theta2_dot)) - 2*r2*cos(theta1 + theta2)*(l1*theta1_dot*sin(theta1) + r2*sin(theta1 + theta2)*(theta1_dot + theta2_dot))))/2 - g*m2*r2*sin(theta1 + theta2);

%display(eq1)
%display(eq2)

EOM = sym('EOM',[2,1]);
EOM(1) = eq1;
EOM(2) = eq2;


EOM_eqlb_points = subs(EOM,[theta1_dot, theta2_dot, theta1_ddot, theta2_ddot, u1, u2], [0,0,0,0,0,0]);

sol = solve(EOM_eqlb_points == 0, [theta1, theta2]);

fprintf("Equilibrium Points");
display(sol.theta1);
display(sol.theta2);


%% Q2.2 (b) Jacobian linearization

m1=1; m2=1; l1=1;l2=1; r1=0.45; r2=0.45; I1=0.084; I2=0.084; g=9.81;

x= sym('x',[4,1]);
x(1) = X1;
x(2) = X2;
x(3) = X3;
x(4) = X4;


x1_dot = X3;
x2_dot = X4;
x3_dot = (I2*u1 - I2*u2 + m2*r2^2*u1*cos(theta1 + theta2)^2 - m2*r2^2*u2*cos(theta1 + theta2)^2 + m2*r2^2*u1*sin(theta1 + theta2)^2 - m2*r2^2*u2*sin(theta1 + theta2)^2 + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2^2*r2^3*theta1_dot^2*cos(theta1 + theta2)^3*sin(theta1) + l1*m2^2*r2^3*theta1_dot^2*sin(theta1 + theta2)^3*cos(theta1) - l1*m2^2*r2^3*theta2_dot^2*cos(theta1 + theta2)^3*sin(theta1) + l1*m2^2*r2^3*theta2_dot^2*sin(theta1 + theta2)^3*cos(theta1) - l1*m2*r2*u2*cos(theta1 + theta2)*cos(theta1) - l1*m2*r2*u2*sin(theta1 + theta2)*sin(theta1) + g*l1*m2^2*r2^2*cos(theta1 + theta2)^2*sin(theta1) - g*l1*m2^2*r2^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)^2 - 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*cos(theta1 + theta2)^3*sin(theta1) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta1 + theta2)^3*cos(theta1) - l1^2*m2^2*r2^2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)*sin(theta1)^2 - I2*l1*m2*r2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1) + I2*l1*m2*r2*theta1_dot^2*sin(theta1 + theta2)*cos(theta1) - I2*l1*m2*r2*theta2_dot^2*cos(theta1 + theta2)*sin(theta1) + I2*l1*m2*r2*theta2_dot^2*sin(theta1 + theta2)*cos(theta1) - l1^2*m2^2*r2^2*theta1_dot^2*cos(theta1 + theta2)^2*cos(theta1)*sin(theta1) + l1^2*m2^2*r2^2*theta1_dot^2*sin(theta1 + theta2)^2*cos(theta1)*sin(theta1) + g*m1*m2*r1*r2^2*cos(theta1 + theta2)^2*sin(theta1) + l1*m2^2*r2^3*theta1_dot^2*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) + l1*m2^2*r2^3*theta2_dot^2*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) + g*m1*m2*r1*r2^2*sin(theta1 + theta2)^2*sin(theta1) - l1*m2^2*r2^3*theta1_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) - l1*m2^2*r2^3*theta2_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) - 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) - 2*I2*l1*m2*r2*theta1_dot*theta2_dot*cos(theta1 + theta2)*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta1 + theta2)*cos(theta1) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1))/(I1*I2 + I1*m2*r2^2*sin(theta1 + theta2)^2 + I2*l1^2*m2*cos(theta1)^2 + I2*m1*r1^2*cos(theta1)^2 + I2*l1^2*m2*sin(theta1)^2 + I2*m1*r1^2*sin(theta1)^2 + I1*m2*r2^2*cos(theta1 + theta2)^2 + l1^2*m2^2*r2^2*cos(theta1 + theta2)^2*sin(theta1)^2 + l1^2*m2^2*r2^2*sin(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*cos(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*cos(theta1 + theta2)^2*sin(theta1)^2 + m1*m2*r1^2*r2^2*sin(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*sin(theta1 + theta2)^2*sin(theta1)^2 - 2*l1^2*m2^2*r2^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)*sin(theta1));
x4_dot = (I1*u2 - I2*u1 + I2*u2 - m2*r2^2*u1*cos(theta1 + theta2)^2 + m2*r2^2*u2*cos(theta1 + theta2)^2 - m2*r2^2*u1*sin(theta1 + theta2)^2 + m2*r2^2*u2*sin(theta1 + theta2)^2 + l1^2*m2*u2*cos(theta1)^2 + m1*r1^2*u2*cos(theta1)^2 + l1^2*m2*u2*sin(theta1)^2 + m1*r1^2*u2*sin(theta1)^2 + I1*g*m2*r2*sin(theta1 + theta2) - I2*g*l1*m2*sin(theta1) - I2*g*m1*r1*sin(theta1) + l1*m2^2*r2^3*theta1_dot^2*cos(theta1 + theta2)^3*sin(theta1) - l1*m2^2*r2^3*theta1_dot^2*sin(theta1 + theta2)^3*cos(theta1) + l1^3*m2^2*r2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1)^3 - l1^3*m2^2*r2*theta1_dot^2*sin(theta1 + theta2)*cos(theta1)^3 + l1*m2^2*r2^3*theta2_dot^2*cos(theta1 + theta2)^3*sin(theta1) - l1*m2^2*r2^3*theta2_dot^2*sin(theta1 + theta2)^3*cos(theta1) - l1*m2*r2*u1*cos(theta1 + theta2)*cos(theta1) + 2*l1*m2*r2*u2*cos(theta1 + theta2)*cos(theta1) - l1*m2*r2*u1*sin(theta1 + theta2)*sin(theta1) + 2*l1*m2*r2*u2*sin(theta1 + theta2)*sin(theta1) - g*l1*m2^2*r2^2*cos(theta1 + theta2)^2*sin(theta1) + g*l1^2*m2^2*r2*sin(theta1 + theta2)*cos(theta1)^2 + l1^3*m2^2*r2*theta1_dot^2*cos(theta1 + theta2)*cos(theta1)^2*sin(theta1) + g*l1*m2^2*r2^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1) - l1^3*m2^2*r2*theta1_dot^2*sin(theta1 + theta2)*cos(theta1)*sin(theta1)^2 - 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)^2 - l1^2*m2^2*r2^2*theta2_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)^2 - g*l1^2*m2^2*r2*cos(theta1 + theta2)*cos(theta1)*sin(theta1) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*cos(theta1 + theta2)^3*sin(theta1) - 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta1 + theta2)^3*cos(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)*sin(theta1)^2 + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)*sin(theta1)^2 + I1*l1*m2*r2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1) - I1*l1*m2*r2*theta1_dot^2*sin(theta1 + theta2)*cos(theta1) + I2*l1*m2*r2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1) - I2*l1*m2*r2*theta1_dot^2*sin(theta1 + theta2)*cos(theta1) + I2*l1*m2*r2*theta2_dot^2*cos(theta1 + theta2)*sin(theta1) - I2*l1*m2*r2*theta2_dot^2*sin(theta1 + theta2)*cos(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta1 + theta2)^2*cos(theta1)*sin(theta1) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta1 + theta2)^2*cos(theta1)*sin(theta1) - 2*l1^2*m2^2*r2^2*theta1_dot^2*sin(theta1 + theta2)^2*cos(theta1)*sin(theta1) - l1^2*m2^2*r2^2*theta2_dot^2*sin(theta1 + theta2)^2*cos(theta1)*sin(theta1) - g*m1*m2*r1*r2^2*cos(theta1 + theta2)^2*sin(theta1) + g*m1*m2*r1^2*r2*sin(theta1 + theta2)*cos(theta1)^2 - l1*m2^2*r2^3*theta1_dot^2*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) - l1*m2^2*r2^3*theta2_dot^2*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) - g*m1*m2*r1*r2^2*sin(theta1 + theta2)^2*sin(theta1) + g*m1*m2*r1^2*r2*sin(theta1 + theta2)*sin(theta1)^2 + l1*m2^2*r2^3*theta1_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) + l1*m2^2*r2^3*theta2_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) - 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)^2 + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta1 + theta2)*sin(theta1 + theta2)*sin(theta1)^2 + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*cos(theta1 + theta2)*sin(theta1) - 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta1 + theta2)*cos(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta1 + theta2)^2*cos(theta1)*sin(theta1) - 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*sin(theta1 + theta2)^2*cos(theta1)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1)^3 - l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta1 + theta2)*cos(theta1)^3 - 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) - g*l1*m1*m2*r1*r2*sin(theta1 + theta2)*sin(theta1)^2 - g*l1*m1*m2*r1*r2*cos(theta1 + theta2)*cos(theta1)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot^2*cos(theta1 + theta2)*cos(theta1)^2*sin(theta1) - l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta1 + theta2)*cos(theta1)*sin(theta1)^2)/(I1*I2 + I1*m2*r2^2*sin(theta1 + theta2)^2 + I2*l1^2*m2*cos(theta1)^2 + I2*m1*r1^2*cos(theta1)^2 + I2*l1^2*m2*sin(theta1)^2 + I2*m1*r1^2*sin(theta1)^2 + I1*m2*r2^2*cos(theta1 + theta2)^2 + l1^2*m2^2*r2^2*cos(theta1 + theta2)^2*sin(theta1)^2 + l1^2*m2^2*r2^2*sin(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*cos(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*cos(theta1 + theta2)^2*sin(theta1)^2 + m1*m2*r1^2*r2^2*sin(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*sin(theta1 + theta2)^2*sin(theta1)^2 - 2*l1^2*m2^2*r2^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)*sin(theta1));

x3_dot = subs(x3_dot,[theta1, theta2, theta1_dot, theta2_dot], [X1,X2,X3,X4]);
x4_dot = subs(x4_dot,[theta1, theta2, theta1_dot, theta2_dot], [X1,X2,X3,X4]);


x_dot= sym('x_dot',[4,1]);
x_dot(1) = x1_dot;
x_dot(2) = x2_dot;
x_dot(3) = x3_dot;
x_dot(4) = x4_dot;

A = jacobian(x_dot,x);
B = jacobian(x_dot,u);

fprintf("A state matrix =");
display(A)
fprintf("B Input matrix =");
display(B)

%% Q2.2 (c) stability properties

fprintf('Linearization around equillibrium point (0,0)');

A = subs(A,[X1,X2,X3,X4], [0,0,0,0]);
B = subs(B,[X1,X2,X3,X4], [0,0,0,0]);

fprintf("A state matrix =");

display(double(A))

fprintf("B input matrix =");
display(double(B))

fprintf('Linearization around equillibrium point (pi,0)');

A1 = subs(A,[X1,X2,X3,X4], [pi,0,0,0]);
B1 = subs(B,[X1,X2,X3,X4], [pi,0,0,0]);

fprintf("A1 state matrix =");

display(simplify(A1))

fprintf("B1 input matrix =");

display(simplify(B1))


fprintf('Linearization around equillibrium point (0,pi)');

A2 = subs(A,[X1,X2,X3,X4], [0,pi,0,0]);
B2 = subs(B,[X1,X2,X3,X4], [0,pi,0,0]);

fprintf("A2 state matrix =");

display(double(A2))

fprintf("B2 input matrix =");

display(double(B2))



%% Q2.2 (d) controllability corresponding to upward position

fprintf("Rank of matrix C =");
rank = rank(ctrb(A,B))

fprintf("Full Rank Matrix, Hence controllable ");
%% Q2.2 (e) state-feedback control design for eigenvalues

eigen_A = eig(double(A))

eigen_A1 = eig(double(A1))

eigen_A2 = eig(double(A2))



fprintf('Chosen values of lambda:')


lambda = [-1,-2,-4,-3]

A = double(A);
B = double(B);


K = place(A,B,lambda)


%% Q2.2 (f) feedback control law function

[t,y]= ode45(@ode_2link,[0 10],[deg2rad(30),deg2rad(45),0,0]);

u=-K*y';


figure('Name','State Trajectories');
subplot(2,2,1)
plot(t,(y(:,1)));
xlabel('t in sec');
ylabel('theta1');
grid on;
subplot(2,2,2)
plot(t,(y(:,2)));
xlabel('t in sec');
ylabel('theta2');
grid on;
subplot(2,2,3)
plot(t,(y(:,3)));
xlabel('t in sec');
ylabel('theta1_dot');
grid on;
subplot(2,2,4)
plot(t,(y(:,4)));
xlabel('t in sec');
ylabel('theta2_dot');
grid on;


figure('Name','Control Inputs Trajectories');
subplot(2,1,1)
plot(t,u(1,:))
xlabel('t in sec');
ylabel('u1');
grid on;

subplot(2,1,2)
plot(t,u(2,:));
xlabel('t in sec');
ylabel('u2');
grid on;
%% Q2.2 (g)
