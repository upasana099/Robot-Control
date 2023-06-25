function dX = ode_planner_2link(t,X)

m1 = 1; m2 = 1; g=9.81; l1 = 1; l2 = 1; r1 = 0.45; r2 = 0.45; 
I1 = 0.084; I2 = 0.084;

dX = zeros(4,1);
X = num2cell(X);
[theta_1, theta_2,q_dot(1), q_dot(2)] = deal(X{:});

% theta_1 = deg2rad(theta_1);
% theta_2 = deg2(theta_2);

tau1 = 0; tau2 = 0;

dX(1) = q_dot(1);
dX(2) = q_dot(2);
dX(3) = (I2*tau1 - I2*tau2 + m2*r2^2*tau1 - m2*r2^2*tau2 + g*l1*m2^2*r2^2*sin(theta_1(t)) + I2*g*l1*m2*sin(theta_1(t)) + I2*g*m1*r1*sin(theta_1(t)) - l1*m2*r2*tau2*cos(theta_2(t)) + l1*m2^2*r2^3*sin(theta_2(t))*diff(theta_1(t), t)^2 + l1*m2^2*r2^3*sin(theta_2(t))*diff(theta_2(t), t)^2 + g*m1*m2*r1*r2^2*sin(theta_1(t)) + 2*l1*m2^2*r2^3*sin(theta_2(t))*diff(theta_1(t), t)*diff(theta_2(t), t) + I2*l1*m2*r2*sin(theta_2(t))*diff(theta_1(t), t)^2 + I2*l1*m2*r2*sin(theta_2(t))*diff(theta_2(t), t)^2 + l1^2*m2^2*r2^2*cos(theta_2(t))*sin(theta_2(t))*diff(theta_1(t), t)^2 - g*l1*m2^2*r2^2*cos(theta_2(t))*sin(theta_1(t) + theta_2(t)) + 2*I2*l1*m2*r2*sin(theta_2(t))*diff(theta_1(t), t)*diff(theta_2(t), t))/(- l1^2*m2^2*r2^2*cos(theta_2(t))^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);
dX(4) = -(I2*tau1 - I1*tau2 - I2*tau2 - l1^2*m2*tau2 - m1*r1^2*tau2 + m2*r2^2*tau1 - m2*r2^2*tau2 + g*l1*m2^2*r2^2*sin(theta_1(t)) + I2*g*l1*m2*sin(theta_1(t)) + I2*g*m1*r1*sin(theta_1(t)) + l1*m2*r2*tau1*cos(theta_2(t)) - 2*l1*m2*r2*tau2*cos(theta_2(t)) - g*l1^2*m2^2*r2*sin(theta_1(t) + theta_2(t)) - I1*g*m2*r2*sin(theta_1(t) + theta_2(t)) + l1*m2^2*r2^3*sin(theta_2(t))*diff(theta_1(t), t)^2 + l1^3*m2^2*r2*sin(theta_2(t))*diff(theta_1(t), t)^2 + l1*m2^2*r2^3*sin(theta_2(t))*diff(theta_2(t), t)^2 + g*m1*m2*r1*r2^2*sin(theta_1(t)) + 2*l1*m2^2*r2^3*sin(theta_2(t))*diff(theta_1(t), t)*diff(theta_2(t), t) + I1*l1*m2*r2*sin(theta_2(t))*diff(theta_1(t), t)^2 + I2*l1*m2*r2*sin(theta_2(t))*diff(theta_1(t), t)^2 + I2*l1*m2*r2*sin(theta_2(t))*diff(theta_2(t), t)^2 + g*l1^2*m2^2*r2*cos(theta_2(t))*sin(theta_1(t)) - g*m1*m2*r1^2*r2*sin(theta_1(t) + theta_2(t)) + 2*l1^2*m2^2*r2^2*cos(theta_2(t))*sin(theta_2(t))*diff(theta_1(t), t)^2 + l1^2*m2^2*r2^2*cos(theta_2(t))*sin(theta_2(t))*diff(theta_2(t), t)^2 - g*l1*m2^2*r2^2*cos(theta_2(t))*sin(theta_1(t) + theta_2(t)) + l1*m1*m2*r1^2*r2*sin(theta_2(t))*diff(theta_1(t), t)^2 + 2*I2*l1*m2*r2*sin(theta_2(t))*diff(theta_1(t), t)*diff(theta_2(t), t) + 2*l1^2*m2^2*r2^2*cos(theta_2(t))*sin(theta_2(t))*diff(theta_1(t), t)*diff(theta_2(t), t) + g*l1*m1*m2*r1*r2*cos(theta_2(t))*sin(theta_1(t)))/(- l1^2*m2^2*r2^2*cos(theta_2(t))^2 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + m1*m2*r1^2*r2^2 + I1*m2*r2^2 + I2*m1*r1^2 + I1*I2);

end