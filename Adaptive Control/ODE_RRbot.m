function dX = ode_rrbot(t,X)
Adaptive_Control = false;

M1 = 1; 
M2 = 1; 
L1 = 1;
L2 = 1;
r1 = 0.45;
r2 = 0.45; 
I1 = 0.084; 
I2 = 0.084; 
g = 9.81; 


dX = zeros(9,1);
X = num2cell(X);
[q1, q2, q1_dot, q2_dot, a1,a2,a3,a4,a5] = deal(X{:});

q1_desired = (pi*t^3)/500 - (3*pi*t^2)/100 - t/18014398509481984 + pi;
q2_desired = (pi*t^3)/1000 - (3*pi*t^2)/200 - t/36028797018963968 + pi/2;

q1_dot_desired = (3*pi*t^2)/500 - (3*pi*t)/50 - 1/18014398509481984;
q2_dot_desired = (3*pi*t^2)/1000 - (3*pi*t)/100 - 1/36028797018963968;

q1_ddot_desired = (3*pi*t)/250 - (3*pi)/50;
q2_ddot_desired = (3*pi*t)/500 - (3*pi)/100;

feed_foward_input = [q1_ddot_desired; q2_ddot_desired];
e = [q1 - q1_desired; q2 - q2_desired];
e_dot = [q1_dot - q1_dot_desired; q2_dot - q2_dot_desired];
% Robust Control input design

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];
l = [-2 -2 -3 -3];

K = place(A,B,l);
Kp = K(:,1:2);            
Kd = K(:,3:4);
O = [0 0; 0 0];
Acl = [O eye(2); -Kp, -Kd] ;
Q = eye(4).*(3.0);

if Adaptive_Control
    P = lyap(Acl',Q);
else
    P = 0;
end
Gamma = eye(5).*25;
G = [e; e_dot];

V = feed_foward_input - Kp*e - Kd*e_dot;

Mmat_hat = [a1 + 2*a2*cos(q2), a3 + a2*cos(q2);a3 + a2*cos(q2),a3];
Cmat_hat = [-a2*q2_dot*sin(q2)*(2*q1_dot + q2_dot);a2*q1_dot^2*sin(q2)];
Gmat_hat = [- a4*g*sin(q1) - a5*g*sin(q1 + q2); -a5*g*sin(q1 + q2)];

U = Mmat_hat.*V + Cmat_hat + Gmat_hat;

tau1 = U(1);
tau2 = U(2);


dX(1) = q1_dot;
dX(2) = q2_dot;
dX(3) = (I2*tau1 - I2*tau2 + M2*r2^2*tau1 - M2*r2^2*tau2 + L1*M2^2*g*r2^2*sin(q1) + I2*L1*M2*g*sin(q1) + I2*M1*g*r1*sin(q1) - L1*M2*r2*tau2*cos(q2) + L1*M2^2*r2^3*q1_dot^2*sin(q2) + L1*M2^2*r2^3*q2_dot^2*sin(q2) + L1^2*M2^2*r2^2*q1_dot^2*cos(q2)*sin(q2) - L1*M2^2*g*r2^2*sin(q1 + q2)*cos(q2) + I2*L1*M2*r2*q1_dot^2*sin(q2) + I2*L1*M2*r2*q2_dot^2*sin(q2) + M1*M2*g*r1*r2^2*sin(q1) + 2*L1*M2^2*r2^3*q1_dot*q2_dot*sin(q2) + 2*I2*L1*M2*r2*q1_dot*q2_dot*sin(q2))/(- L1^2*M2^2*r2^2*cos(q2)^2 + L1^2*M2^2*r2^2 + I2*L1^2*M2 + M1*M2*r1^2*r2^2 + I1*M2*r2^2 + I2*M1*r1^2 + I1*I2);
dX(4) = -(I2*tau1 - I1*tau2 - I2*tau2 - L1^2*M2*tau2 - M1*r1^2*tau2 + M2*r2^2*tau1 - M2*r2^2*tau2 - L1^2*M2^2*g*r2*sin(q1 + q2) + L1*M2^2*g*r2^2*sin(q1) - I1*M2*g*r2*sin(q1 + q2) + I2*L1*M2*g*sin(q1) + I2*M1*g*r1*sin(q1) + L1*M2*r2*tau1*cos(q2) - 2*L1*M2*r2*tau2*cos(q2) + L1*M2^2*r2^3*q1_dot^2*sin(q2) + L1^3*M2^2*r2*q1_dot^2*sin(q2) + L1*M2^2*r2^3*q2_dot^2*sin(q2) + 2*L1^2*M2^2*r2^2*q1_dot^2*cos(q2)*sin(q2) + L1^2*M2^2*r2^2*q2_dot^2*cos(q2)*sin(q2) - L1*M2^2*g*r2^2*sin(q1 + q2)*cos(q2) + L1^2*M2^2*g*r2*cos(q2)*sin(q1) - M1*M2*g*r1^2*r2*sin(q1 + q2) + I1*L1*M2*r2*q1_dot^2*sin(q2) + I2*L1*M2*r2*q1_dot^2*sin(q2) + I2*L1*M2*r2*q2_dot^2*sin(q2) + M1*M2*g*r1*r2^2*sin(q1) + 2*L1*M2^2*r2^3*q1_dot*q2_dot*sin(q2) + 2*L1^2*M2^2*r2^2*q1_dot*q2_dot*cos(q2)*sin(q2) + L1*M1*M2*r1^2*r2*q1_dot^2*sin(q2) + 2*I2*L1*M2*r2*q1_dot*q2_dot*sin(q2) + L1*M1*M2*g*r1*r2*cos(q2)*sin(q1))/(- L1^2*M2^2*r2^2*cos(q2)^2 + L1^2*M2^2*r2^2 + I2*L1^2*M2 + M1*M2*r1^2*r2^2 + I1*M2*r2^2 + I2*M1*r1^2 + I1*I2);

q1_ddot = dX(3);
q2_ddot = dX(4);

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

Phi = Mmat_hat\Y;
alpha_tilda = -(Gamma\(Phi'*B'*P*G));

dX(5) = alpha_tilda(1);
dX(6) = alpha_tilda(2);
dX(7) = alpha_tilda(3);
dX(8) = alpha_tilda(4);
dX(9) = alpha_tilda(5);
end
