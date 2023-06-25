set(0, 'defaultFigureRenderer', 'painters')
clear; close; clc;
% ROS Setup
rosinit;

j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');
JointStates = rossubscriber('/rrbot/joint_states');

tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(200), deg2rad(125)];
resp = call(client,req,'Timeout',3);

tic;
t = 0;
i = 1;
j=1;

while(t < 10)
    t = toc;
    time(i) = t;
    % read the joint states
    
    jointData = receive(JointStates);
    
    % inspect the "jointData" variable in MATLAB to get familiar with its structure
    
    % implement your state feedback controller below

    X=[jointData.Position(1);jointData.Position(2);jointData.Velocity(1);jointData.Velocity(2)];
    
    theta_1(i) = jointData.Position(1);
    theta_1_dot(i) = jointData.Velocity(1);
    theta_2(i) = jointData.Position(2);
    theta_2_dot(i) = jointData.Velocity(2);
    
 
    while (theta_1(i)>2*pi)
        theta_1(i)=theta_1(i)-2*pi;
    end
    while (theta_2(i)>2*pi)
        theta_2(i)=theta_2(i)-2*pi;
    end
    while (theta_1(i)<-0.01)
        theta_1(i)=theta_1(i)+2*pi;
    end
    while (theta_2(i)<-0.01)
        theta_2(i)=theta_2(i)+2*pi;
    end


    X=[theta_1(i);theta_2(i);theta_1_dot(i);theta_2_dot(i)];

    [~,u] = robustcon(t,X);

    tau1.Data = u(1);
    tau2.Data = u(2);

    u1(i) = tau1.Data;
    u2(i) = tau2.Data;
    send(j1_effort,tau1);
    send(j2_effort,tau2);
    
    i = i+1;
    % sample the time, joint state values, and calculated torques here to be plotted at the end   
   


    qd_arr(j,1)=(pi*t^3)/500 - (3*pi*t^2)/100 + pi;
    qd_arr(j,2)=(pi*t^3)/1000 - (3*pi*t^2)/200 + pi/2;
    j=j+1;

end

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

% disconnect from roscore
rosshutdown;

% Plot the trajectories


figure;
plot(time,u1)
xlabel('time in sec');
ylabel('u1');


figure;
plot(time,u2)
xlabel('time in sec');    
ylabel('u2');

figure;
%hold on;
plot(time,theta_1);
plot(time, qd_arr(:,1));
xlabel('time in sec');
ylabel('theta 1');
%hold off;

figure;

plot(time,theta_1_dot);
xlabel('time t in sec');
ylabel('theta1 dot');


figure;
%hold on;
plot(time,theta_2);
plot(time, qd_arr(:,2));
xlabel('time t in sec');
ylabel('theta 2');
%hold off;


figure;
%hold on;
plot(time,theta_2_dot);
xlabel('time t in sec');
ylabel('theta2 dot');
%hold off;
