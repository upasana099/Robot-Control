%% RBE 502 : Prog Assignment 4
% Upasana Mahanti

clear; close; clc;
syms theta1 theta2 theta1_dot theta2_dot theta1_Ddot theta2_Ddot m1 m2 l1 l2 g I1 I2 r1 r2 u1 u2 m1_hat m2_hat I1_hat I2_hat 'real'; 
syms theta1_d theta2_d theta1_dot_d theta2_dot_d q q1 q2 dq1 dq2 
syms a0 a1 a2 a3 b0 b1 b2 b3 t x1 x2 x3 x4 

%% Part (a) : Cubic Polynomial Trajectory
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

% disp([A1';B1'])   
qd=[A1';B1']*[1;t;t^2;t^3];
qd_dot=[A1';B1']*[0;1;2*t;3*t^2];
qd_ddot=[A1';B1']*[0;0;2;6*t];

% 
% figure
% subplot(2,3,1)
% fplot(qd(1), [0 10])
% xlabel('time t in sec');
% ylabel('theta1.d');
% subplot(2,3,2)
% fplot(qd_dot(1), [0 10])
% xlabel('time t in sec');
% ylabel('theta1.d.dot');
% subplot(2,3,3)
% fplot(qd_ddot(1), [0 10])
% xlabel('time t in sec');
% ylabel('theta1.d.ddot');
% subplot(2,3,4)
% fplot(qd(2), [0 10])
% xlabel('time t in sec');
% ylabel('theta2.d');
% subplot(2,3,5)
% fplot(qd_dot(2), [0 10])
% xlabel('time t in sec');
% ylabel('theta2.d.dot');
% subplot(2,3,6)
% fplot(qd_ddot(2), [0 10])
% xlabel('time t in sec');
% ylabel('theta2.d.ddot');


%% Part (b) : Manipulator Form

u = [u1;u2];


a = I1 + I2 + m1*r1^2 + m2*(l1^2 + r2^2);
b = m2*l1*r2;
d = I2 + m2*r2^2;

Mmat= [a+2*b*cos(q2), d+b*cos(q2); d+b*cos(q2), d];
Cmat= [-b*sin(q2)*dq2, -b*sin(q2)*(dq1+dq2); b*sin(q2)*dq1,0];
Gmat= [-m1*g*r1*sin(q1)-m2*g*(l1*sin(q1)+r2*sin(q1+q2)); -m2*g*r2*sin(q1+q2)];


Mmat=subs(Mmat,[q1,q2,dq1,dq2],[theta1, theta2, theta1_dot, theta2_dot]);
Cmat=subs(Cmat,[q1,q2,dq1,dq2],[theta1, theta2, theta1_dot, theta2_dot]);
Gmat=subs(Gmat,[q1,q2,dq1,dq2],[theta1, theta2, theta1_dot, theta2_dot]);


Mmat_hat=subs(Mmat,[m1,m2,I1,I2],[m1_hat,m2_hat,I1_hat,I2_hat]);
Cmat_hat=subs(Cmat,[m1,m2,I1,I2],[m1_hat,m2_hat,I1_hat,I2_hat]);
Gmat_hat=subs(Gmat,[m1,m2,I1,I2],[m1_hat,m2_hat,I1_hat,I2_hat]);


% 

X=sym('X',[4,1]);
X(1) = theta1;
X(2) = theta2;
X(3) = theta1_dot; 
X(4) = theta2_dot;

Xd=[qd;qd_dot];
e=X-Xd;

%% Part (c): Robust Inverse dynamics control law
% Eigenvalue placement method

A=[0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];
B=[0,0;0,0;1,0;0,1];
p=[-1,-1,-2,-2];
K=place(A,B,p);

% Lyapunov matrix

Acl= A-B*K;
%Q=eye(4).*10;
Q = diag([30,30,10,10])

P=lyap(Acl',Q)


% 

[t,y] = ode45(@robustcon,[0 10],[deg2rad(200),deg2rad(125),0,0]); % angles in radian

u=zeros([length(t),2]);

for i=1:length(t)
    [~,u(i,:)] = robustcon(t(i,1),y(i,:));
end


% Plots

figure;

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
ylabel('theta1.dot');
hold off;

subplot(2,2,4)
hold on;
fplot(qd_dot(2), [0 10])
plot(t,y(:,4));
xlabel('time t in sec');
ylabel('theta2.dot');
hold off;


figure;
subplot(2,1,1)
plot(t,u(:,1))
xlabel('time t in sec');
ylabel('u1');

subplot(2,1,2)
plot(t,u(:,2));
xlabel('time t in sec');
ylabel('u2');

