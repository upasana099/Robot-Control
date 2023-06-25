% Upasana Mahanti
clc;
clear all;
close all;
% System Parameters

M1 = 1; 
M2 = 1; 
L1 = 1; 
L2 = 1; 
r1 = 0.45; 
r2 = 0.45; 
I1 = 0.084; 
I2 = 0.084; 
g = 9.81; 

% M1_hat = 0.75; 
% M2_hat = 0.75; 
% I1_hat = 0.063; 
% I2_hat = 0.063;

syms q1 q2 q1_dot q2_dot q1_ddot q2_ddot tau1 tau2 g 'real'
syms M1 M2 M1_hat M2_hat I1 I2 I1_hat I2_hat L1 L2 r1 r2 'real'
syms x1_dot y1_dot x2_dot y2_dot 'real'
syms a1 a2 a3 a4 a5 'real'
syms t 'real'

% State Space Representation
X = sym('X', [4,1]);
X(1) = q1;
X(2) = q2;
X(3) = q1_dot;
X(4) = q2_dot;

% Dynamic Equations

x1_dot = (q1_dot)*(r1)*(cos(q1));
y1_dot = -(q1_dot)*(r1)*(sin(q1));

x2_dot = (q1_dot)*(L1)*(cos(q1)) + (q1_dot + q2_dot)*(r2)*(cos(q1 + q2));
y2_dot = -(q1_dot)*(L1)*(sin(q1)) - (q1_dot + q2_dot)*(r2)*(sin(q1 + q2));

KE1 = (1/2)*(I1)*(q1_dot*q1_dot) + (1/2)*(M1)*((x1_dot*x1_dot) + (y1_dot*y1_dot));
KE2 = (1/2)*(I2)*((q2_dot + q1_dot)*(q2_dot + q1_dot)) + (1/2)*(M2)*((x2_dot*x2_dot) + (y2_dot*y2_dot));

PE1 = M1*g*r1*cos(q1);
PE2 = M2*g*(L1*(cos(q1)) + r2*(cos(q1 + q2)));
L = KE1 + KE2 - PE1 - PE2;

u = [tau1;tau2];
q = [q1;q2];
dq = [q1_dot; q2_dot];
ddq = [q1_ddot; q2_ddot];

DL_Dq = gradient(L,q);  
DL_Ddq = gradient(L,dq); 
dDL_dtDdq = jacobian(DL_Ddq,[q;dq])*[dq;ddq];

EOM = simplify(dDL_dtDdq - DL_Dq -u);

EOM_numerical = subs(EOM,[M1,M2,L1,L2,r1,r2,I1,I2,g],[1,1,1,1,0.45,0.45,0.084,0.084,9.81]);


% Feedback Linearization

a = I1 + I2 + M1*r1^2 + M2*(L1^2 + r2^2);
b = M2*L1*r2;
d = I2 + M2*r2^2;

Mmat= [a+2*b*cos(q2), d+b*cos(q2); d+b*cos(q2), d];
Cmat= [-b*sin(q2)*q2_dot, -b*sin(q2)*(q1_dot + q2_dot); b*sin(q2)*q1_dot,0];
Gmat= [-M1*g*r1*sin(q1)-M2*g*(L1*sin(q1)+r2*sin(q1+q2)); -M2*g*r2*sin(q1+q2)];


EOM_FL = simplify(Mmat*[q1_ddot; q2_ddot] + Cmat*[q1_dot; q2_dot] + Gmat)


% Linear Parameteric Form

Y = [q1_ddot, ...
     cos(q2)*(2*q1_ddot + q2_ddot) - 2*sin(q2)*q1_dot*q2_dot - sin(q2)*q2_dot^2, ...
     q2_ddot, ...
    -sin(q1)*g, ...
    -sin(q1 + q2)*g; ...
     0, ...
     sin(q2)*q1_dot^2 + cos(q2)*q1_ddot, ...
     q1_ddot + q2_ddot, ...
     0, ...
     -sin(q1+q2)*g];


alpha = [M2*L1^2 + M1*r1^2 + M2*r2^2 + I1 + I2;
            M2*L1*r2;
            M2*r2^2 + I2;
            M1*r1 + M2*L1;
            M2*r2];


tau_Yalpha = Y * alpha;


% tau : manipulator equation
tau = Mmat * [q1_ddot; q2_ddot] + Cmat * [q1_dot; q2_dot] + Gmat;

% Compare the terms in the product Y * alpha with the manipulator form equations
result1 = simplify(tau_Yalpha(1) - tau(1));
result2 = simplify(tau_Yalpha(2) - tau(2));

% checking to see if Y matrix and alpha vector match the manipulator form equations
if result1 == 0 && result2 == 0
    disp('The provided Y matrix and alpha vector match the manipulator form equations.')
else
    disp('The provided Y matrix and alpha vector do not match the manipulator form equations.')
end


EOM_LP = simplify(Y*alpha)

alpha_para = [a1;a2;a3;a4;a5];

EOM_LP_1 = Y*alpha_para

g_q_symbolic = simplify(subs(EOM_LP_1,[q1_dot, q2_dot, q1_ddot, q2_ddot, tau1, tau2],[0,0,0,0,0,0]))


M_q_temp_s = subs(EOM_LP_1,[q1_dot, q2_dot, tau1, tau2],[0,0,0,0]) - g_q_symbolic;
M_q_symbolic = simplify(jacobian(M_q_temp_s,[q1_ddot,q2_ddot]))


C_q_qd_symbolic = simplify(subs(EOM_LP_1 - g_q_symbolic - M_q_temp_s,[tau1, tau2],[0,0]))


