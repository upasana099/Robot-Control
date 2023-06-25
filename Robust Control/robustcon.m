function [dx,u] = robustcon(t,x)

m1=1; m2=1; l1=1; l2=1; r1=0.45; r2=0.45; I1=0.084; I2=0.084; g=9.81;
m1_hat=0.75; m2_hat=0.75; I1_hat=0.063; I2_hat=0.063;

dx = zeros(4,1);
x = num2cell(x);
[theta1, theta2, theta1_dot, theta2_dot] = deal(x{:});


e = [- (pi*t^3)/500 + (3*pi*t^2)/100 - pi + theta1; -(pi*t^3)/1000 + (3*pi*t^2)/200 - pi/2 + theta2; theta1_dot + (3*pi*t)/50 - (3*pi*t^2)/500; theta2_dot + (3*pi*t)/100 - (3*pi*t^2)/1000];
qd_ddot=[(3*pi*t)/250 - (3*pi)/50;(3*pi*t)/500 - (3*pi)/100];

A=[0,0,1,0;0,0,0,1;0,0,0,0;0,0,0,0];
B=[0,0;0,0;1,0;0,1];

K = [2.0000 , 0 ,     3.0000, 0; 
     0,       2.0000, 0,      3.0000];


P =[95.0000         0   25.0000         0;
         0   95.0000         0   25.0000;
   25.0000         0   10.0000         0;
         0   25.0000         0   10.0000];




Mmat_hat=[I1_hat + I2_hat + m1_hat*r1^2 + m2_hat*(l1^2 + r2^2) + 2*l1*m2_hat*r2*cos(theta2), m2_hat*r2^2 + l1*m2_hat*cos(theta2)*r2 + I2_hat; m2_hat*r2^2 + l1*m2_hat*cos(theta2)*r2 + I2_hat,m2_hat*r2^2 + I2_hat];
Cmat_hat = [-l1*m2_hat*r2*theta2_dot*sin(theta2), -l1*m2_hat*r2*sin(theta2)*(theta1_dot + theta2_dot);l1*m2_hat*r2*theta1_dot*sin(theta2),0];
Gmat_hat=[- g*m2_hat*(r2*sin(theta1 + theta2) + l1*sin(theta1)) - g*m1_hat*r1*sin(theta1);-g*m2_hat*r2*sin(theta1 + theta2)];


% Uncertainty Upper bound
rho = 7.0;


% Robust control term
% vr=0;
                                                                                                                                                                                                                                                                                                                                                                                            
phi=4.0;

if phi > 0
    if norm(B'*P*e) > phi
        vr= -rho*(B'*P*e)/(norm(B'*P*e));
    else
        vr= -rho*(B'*P*e)/(phi);
    end
else
    if norm(B'*P*e) ~= 0
        vr= -rho*(B'*P*e)/(norm(B'*P*e));
    else
        vr= 0;
    end
end
                                                                                                                                    
% Virtual control input v
v = qd_ddot - K*e + vr;

u = Mmat_hat*v + Cmat_hat*[theta1_dot;theta2_dot] + Gmat_hat;


u1=u(1);
u2=u(2);

dx(1) = theta1_dot;
dx(2) = theta2_dot;
dx(3) = (I2*u1 - I2*u2 + m2*r2^2*u1*cos(theta1 + theta2)^2 - m2*r2^2*u2*cos(theta1 + theta2)^2 + m2*r2^2*u1*sin(theta1 + theta2)^2 - m2*r2^2*u2*sin(theta1 + theta2)^2 + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) - l1*m2^2*r2^3*theta1_dot^2*cos(theta1 + theta2)^3*sin(theta1) + l1*m2^2*r2^3*theta1_dot^2*sin(theta1 + theta2)^3*cos(theta1) - l1*m2^2*r2^3*theta2_dot^2*cos(theta1 + theta2)^3*sin(theta1) + l1*m2^2*r2^3*theta2_dot^2*sin(theta1 + theta2)^3*cos(theta1) - l1*m2*r2*u2*cos(theta1 + theta2)*cos(theta1) - l1*m2*r2*u2*sin(theta1 + theta2)*sin(theta1) + g*l1*m2^2*r2^2*cos(theta1 + theta2)^2*sin(theta1) - g*l1*m2^2*r2^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1) + l1^2*m2^2*r2^2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)^2 - 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*cos(theta1 + theta2)^3*sin(theta1) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta1 + theta2)^3*cos(theta1) - l1^2*m2^2*r2^2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)*sin(theta1)^2 - I2*l1*m2*r2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1) + I2*l1*m2*r2*theta1_dot^2*sin(theta1 + theta2)*cos(theta1) - I2*l1*m2*r2*theta2_dot^2*cos(theta1 + theta2)*sin(theta1) + I2*l1*m2*r2*theta2_dot^2*sin(theta1 + theta2)*cos(theta1) - l1^2*m2^2*r2^2*theta1_dot^2*cos(theta1 + theta2)^2*cos(theta1)*sin(theta1) + l1^2*m2^2*r2^2*theta1_dot^2*sin(theta1 + theta2)^2*cos(theta1)*sin(theta1) + g*m1*m2*r1*r2^2*cos(theta1 + theta2)^2*sin(theta1) + l1*m2^2*r2^3*theta1_dot^2*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) + l1*m2^2*r2^3*theta2_dot^2*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) + g*m1*m2*r1*r2^2*sin(theta1 + theta2)^2*sin(theta1) - l1*m2^2*r2^3*theta1_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) - l1*m2^2*r2^3*theta2_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) - 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) - 2*I2*l1*m2*r2*theta1_dot*theta2_dot*cos(theta1 + theta2)*sin(theta1) + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta1 + theta2)*cos(theta1) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1))/(I1*I2 + I1*m2*r2^2*sin(theta1 + theta2)^2 + I2*l1^2*m2*cos(theta1)^2 + I2*m1*r1^2*cos(theta1)^2 + I2*l1^2*m2*sin(theta1)^2 + I2*m1*r1^2*sin(theta1)^2 + I1*m2*r2^2*cos(theta1 + theta2)^2 + l1^2*m2^2*r2^2*cos(theta1 + theta2)^2*sin(theta1)^2 + l1^2*m2^2*r2^2*sin(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*cos(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*cos(theta1 + theta2)^2*sin(theta1)^2 + m1*m2*r1^2*r2^2*sin(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*sin(theta1 + theta2)^2*sin(theta1)^2 - 2*l1^2*m2^2*r2^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)*sin(theta1));
dx(4) = (I1*u2 - I2*u1 + I2*u2 - m2*r2^2*u1*cos(theta1 + theta2)^2 + m2*r2^2*u2*cos(theta1 + theta2)^2 - m2*r2^2*u1*sin(theta1 + theta2)^2 + m2*r2^2*u2*sin(theta1 + theta2)^2 + l1^2*m2*u2*cos(theta1)^2 + m1*r1^2*u2*cos(theta1)^2 + l1^2*m2*u2*sin(theta1)^2 + m1*r1^2*u2*sin(theta1)^2 + I1*g*m2*r2*sin(theta1 + theta2) - I2*g*l1*m2*sin(theta1) - I2*g*m1*r1*sin(theta1) + l1*m2^2*r2^3*theta1_dot^2*cos(theta1 + theta2)^3*sin(theta1) - l1*m2^2*r2^3*theta1_dot^2*sin(theta1 + theta2)^3*cos(theta1) + l1^3*m2^2*r2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1)^3 - l1^3*m2^2*r2*theta1_dot^2*sin(theta1 + theta2)*cos(theta1)^3 + l1*m2^2*r2^3*theta2_dot^2*cos(theta1 + theta2)^3*sin(theta1) - l1*m2^2*r2^3*theta2_dot^2*sin(theta1 + theta2)^3*cos(theta1) - l1*m2*r2*u1*cos(theta1 + theta2)*cos(theta1) + 2*l1*m2*r2*u2*cos(theta1 + theta2)*cos(theta1) - l1*m2*r2*u1*sin(theta1 + theta2)*sin(theta1) + 2*l1*m2*r2*u2*sin(theta1 + theta2)*sin(theta1) - g*l1*m2^2*r2^2*cos(theta1 + theta2)^2*sin(theta1) + g*l1^2*m2^2*r2*sin(theta1 + theta2)*cos(theta1)^2 + l1^3*m2^2*r2*theta1_dot^2*cos(theta1 + theta2)*cos(theta1)^2*sin(theta1) + g*l1*m2^2*r2^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1) - l1^3*m2^2*r2*theta1_dot^2*sin(theta1 + theta2)*cos(theta1)*sin(theta1)^2 - 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)^2 - l1^2*m2^2*r2^2*theta2_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)^2 - g*l1^2*m2^2*r2*cos(theta1 + theta2)*cos(theta1)*sin(theta1) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*cos(theta1 + theta2)^3*sin(theta1) - 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*sin(theta1 + theta2)^3*cos(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)*sin(theta1)^2 + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)*sin(theta1)^2 + I1*l1*m2*r2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1) - I1*l1*m2*r2*theta1_dot^2*sin(theta1 + theta2)*cos(theta1) + I2*l1*m2*r2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1) - I2*l1*m2*r2*theta1_dot^2*sin(theta1 + theta2)*cos(theta1) + I2*l1*m2*r2*theta2_dot^2*cos(theta1 + theta2)*sin(theta1) - I2*l1*m2*r2*theta2_dot^2*sin(theta1 + theta2)*cos(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot^2*cos(theta1 + theta2)^2*cos(theta1)*sin(theta1) + l1^2*m2^2*r2^2*theta2_dot^2*cos(theta1 + theta2)^2*cos(theta1)*sin(theta1) - 2*l1^2*m2^2*r2^2*theta1_dot^2*sin(theta1 + theta2)^2*cos(theta1)*sin(theta1) - l1^2*m2^2*r2^2*theta2_dot^2*sin(theta1 + theta2)^2*cos(theta1)*sin(theta1) - g*m1*m2*r1*r2^2*cos(theta1 + theta2)^2*sin(theta1) + g*m1*m2*r1^2*r2*sin(theta1 + theta2)*cos(theta1)^2 - l1*m2^2*r2^3*theta1_dot^2*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) - l1*m2^2*r2^3*theta2_dot^2*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) - g*m1*m2*r1*r2^2*sin(theta1 + theta2)^2*sin(theta1) + g*m1*m2*r1^2*r2*sin(theta1 + theta2)*sin(theta1)^2 + l1*m2^2*r2^3*theta1_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) + l1*m2^2*r2^3*theta2_dot^2*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) + 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*cos(theta1 + theta2)*sin(theta1 + theta2)^2*sin(theta1) - 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)^2 + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta1 + theta2)*sin(theta1 + theta2)*sin(theta1)^2 + 2*I2*l1*m2*r2*theta1_dot*theta2_dot*cos(theta1 + theta2)*sin(theta1) - 2*I2*l1*m2*r2*theta1_dot*theta2_dot*sin(theta1 + theta2)*cos(theta1) + 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*cos(theta1 + theta2)^2*cos(theta1)*sin(theta1) - 2*l1^2*m2^2*r2^2*theta1_dot*theta2_dot*sin(theta1 + theta2)^2*cos(theta1)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot^2*cos(theta1 + theta2)*sin(theta1)^3 - l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta1 + theta2)*cos(theta1)^3 - 2*l1*m2^2*r2^3*theta1_dot*theta2_dot*cos(theta1 + theta2)^2*sin(theta1 + theta2)*cos(theta1) - g*l1*m1*m2*r1*r2*sin(theta1 + theta2)*sin(theta1)^2 - g*l1*m1*m2*r1*r2*cos(theta1 + theta2)*cos(theta1)*sin(theta1) + l1*m1*m2*r1^2*r2*theta1_dot^2*cos(theta1 + theta2)*cos(theta1)^2*sin(theta1) - l1*m1*m2*r1^2*r2*theta1_dot^2*sin(theta1 + theta2)*cos(theta1)*sin(theta1)^2)/(I1*I2 + I1*m2*r2^2*sin(theta1 + theta2)^2 + I2*l1^2*m2*cos(theta1)^2 + I2*m1*r1^2*cos(theta1)^2 + I2*l1^2*m2*sin(theta1)^2 + I2*m1*r1^2*sin(theta1)^2 + I1*m2*r2^2*cos(theta1 + theta2)^2 + l1^2*m2^2*r2^2*cos(theta1 + theta2)^2*sin(theta1)^2 + l1^2*m2^2*r2^2*sin(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*cos(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*cos(theta1 + theta2)^2*sin(theta1)^2 + m1*m2*r1^2*r2^2*sin(theta1 + theta2)^2*cos(theta1)^2 + m1*m2*r1^2*r2^2*sin(theta1 + theta2)^2*sin(theta1)^2 - 2*l1^2*m2^2*r2^2*cos(theta1 + theta2)*sin(theta1 + theta2)*cos(theta1)*sin(theta1));


end