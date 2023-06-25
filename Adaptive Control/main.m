clc; clear; close;
syms q1 q2 q1_dot q2_dot q1_ddot q2_ddot tau1 tau2 g 'real'
syms t 'real'
syms a1 a2 a3 a4 a5 'real'
visualization = true;
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

aplha_desired = [M2*L1^2 + M1*r1^2 + M2*r2^2 + I1 + I2;
M2*L1*r2;
M2*r2^2 + I2;
M1*r1 + M2*L1;
M2*r2]

% Robust Inverse dynamics control law

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];

l = [-2 -2 -3 -3];

K = place(A,B,l)
Kp = K(:,1:2)             
Kd = K(:,3:4)
O = [0 0; 0 0];
Acl = [O eye(2); -Kp, -Kd]
Q = eye(4).*(3.0)
% P = lyap(Acl',Q)

if Adaptive_Control
    P = lyap(Acl',Q);
else
    P = 0;
end
Gamma = eye(5).*25
T = 10;
q1_initial = 200;
q2_initial = 125;
alpha0 = [M2*L1^2 + M1*r1^2 + M2*r2^2 + I1 + I2;
            M2*L1*r2;
            M2*r2^2 + I2;
            M1*r1 + M2*L1;
            M2*r2].*0.75; 
y0 = [deg2rad(q1_initial), deg2rad(q2_initial) 0, 0,alpha0'];

[time,y] = ode45(@ode_rrbot,[0,T],y0);

q1_desired = (pi*time.^3)/500 - (3*pi*time.^2)/100 - time/18014398509481984 + pi;
q2_desired = (pi*time.^3)/1000 - (3*pi*time.^2)/200 - time/36028797018963968 + pi/2;

q1_dot_desired = (3*pi*time.^2)/500 - (3*pi*time)/50 - 1/18014398509481984;
q2_dot_desired = (3*pi*time.^2)/1000 - (3*pi*time)/100 - 1/36028797018963968;

q1_ddot_desired = (3*pi*time)/250 - (3*pi)/50;
q2_ddot_desired = (3*pi*time)/500 - (3*pi)/100;




for i = 1:size(time)
    q1 = y(i,1);
    q2 = y(i,2);
    q1_dot = y(i,3);
    q2_dot = y(i,4);
    a1 = y(i,5);
    a2 = y(i,6);
    a3 = y(i,7);
    a4 = y(i,8);
    a5 = y(i,9);
    
    q1_desired_1 = (pi*time(i)^3)/500 - (3*pi*time(i)^2)/100 - time(i)/18014398509481984 + pi;
    q2_desired_1 = (pi*time(i)^3)/1000 - (3*pi*time(i)^2)/200 - time(i)/36028797018963968 + pi/2;
    q1_dot_desired_1 = (3*pi*time(i)^2)/500 - (3*pi*time(i))/50 - 1/18014398509481984;
    q2_dot_desired_1 = (3*pi*time(i)^2)/1000 - (3*pi*time(i))/100 - 1/36028797018963968;
    q1_ddot_desired_1 = (3*pi*time(i))/250 - (3*pi)/50;
    q2_ddot_desired_1 = (3*pi*time(i))/500 - (3*pi)/100;

    feed_foward_input = [q1_ddot_desired_1; q2_ddot_desired_1];
    e = [y(i,1) - q1_desired_1; y(i,2) - q2_desired_1];
    e_dot = [y(i,3) - q1_dot_desired_1; y(i,4) - q2_dot_desired_1];

    G = [e; e_dot];
    
    Mmat_hat = [a1 + 2*a2*cos(q2), a3 + a2*cos(q2);a3 + a2*cos(q2),                  a3];
    Cmat_hat = [-a2*q2_dot*sin(q2)*(2*q1_dot + q2_dot);a2*q1_dot^2*sin(q2)];
    Gmat_hat = [- a4*g*sin(q1) - a5*g*sin(q1 + q2); -a5*g*sin(q1 + q2)];

    V = feed_foward_input - Kp*e - Kd*e_dot;
    U = Mmat_hat*V + Cmat_hat + Gmat_hat;
    U1(i) = U(1);
    U2(i) = U(2);
end


% plots
figure
hold on
subplot(2,2,1)
plot(time,y(:,1))
hold on
plot(time, q1_desired,'Color','red','LineStyle','--')
xlabel('time')
ylabel('rad')
title('q1')

subplot(2,2,2)
plot(time,y(:,2))
hold on
plot(time, q2_desired,'Color','red','LineStyle','--')
xlabel('time')
ylabel('rad')
title('q2')

subplot(2,2,3)
plot(time,y(:,3))
hold on
plot(time, q1_dot_desired,'Color','red','LineStyle','--')
xlabel('time')
ylabel('rad/s')
title('q1-dot')

subplot(2,2,4)
plot(time,y(:,4))
hold on
plot(time, q2_dot_desired,'Color','red','LineStyle','--')
xlabel('time')
ylabel('rad/s')
title('q2-dot')
hold off

figure
hold on
subplot(2,1,1)
plot(time,U1)
xlabel('time')
ylabel('Nm')
title('tau1')

subplot(2,1,2)
plot(time,U2)
xlabel('time')
ylabel('Nm')
title('tau2')
hold off

figure
subplot(2,3,1)
plot(time,y(:,5))
hold on
xlabel('time')
ylabel('a1')
title('a1')
hold off

subplot(2,3,2)
plot(time,y(:,6))
hold on
xlabel('time')
ylabel('a2')
title('a2')
hold off

subplot(2,3,3)
plot(time,y(:,7))
hold on
xlabel('time')
ylabel('a3')
title('a3')
hold off

subplot(2,3,4)
plot(time,y(:,8))
hold on
xlabel('time')
ylabel('a4')
title('a4')
hold off

subplot(2,3,5)
plot(time,y(:,9))
hold on
xlabel('time')
ylabel('a5')
title('a5')
hold off
pause(10)
close all
% visualization
figure
if(visualization)

    q1_desired_plot = (pi*time.^3)/500 - (3*pi*time.^2)/100 - time/18014398509481984 + pi;
    q2_desired_plot = (pi*time.^3)/1000 - (3*pi*time.^2)/200 - time/36028797018963968 + pi/2;
   
    l1 = 1;
    l2 = 1;

    x1_desired_plot = l1*sin(q1_desired_plot);
    y1_desired_plot = l1*cos(q1_desired_plot);

    x2_desired_plot = l1*sin(q1_desired_plot) + l2*sin(q1_desired_plot+q2_desired_plot);
    y2_desired_plot = l1*cos(q1_desired_plot) + l2*cos(q1_desired_plot+q2_desired_plot);
    

    x1_pos= l1*sin(y(:,1));
    x2_pos= l1*sin(y(:,1)) + l2*sin(y(:,1)+y(:,2));
    y1_pos= l1*cos(y(:,1));
    y2_pos= l1*cos(y(:,1)) + l2*cos(y(:,1)+y(:,2));

    for i=1:size(y)
        plot(x1_desired_plot, y1_desired_plot,'Color','red','LineStyle','--')
        hold on
        plot(x2_desired_plot, y2_desired_plot,'Color','red','LineStyle','--')
        hold on
        plot([0 x1_pos(i) x2_pos(i)],[0 y1_pos(i) y2_pos(i)],'black')
        hold on
        plot(x2_pos(1:i),y2_pos(1:i),'black')
        hold on
        plot(x1_pos(1:i),y1_pos(1:i),'black')
        hold on
        plot(x1_pos(i),y1_pos(i),'.')
        xlim([-2.25 2.25])
        ylim([-2.25 2.25])
        pause(0.000075)
        hold off
    end
end
