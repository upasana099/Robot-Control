clear; close; clc;
% ROS Setup
rosshutdown;
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
req.JointPositions = [deg2rad(30), deg2rad(45)];
resp = call(client,req,'Timeout',3);
tic;
t = 0;
i = 1;

while(t < 10)
    t = toc;
    time(i) = t;
    % read the joint states
    jointData = receive(JointStates);
    % inspect the "jointData" variable in MATLAB to get familiar with its structure
    % implement your state feedback controller below
    
    K=[25.7954,8.2635,9.7870,3.6006; 6.4934,6.0896,2.8937,1.4718];
    X=[jointData.Position(1);jointData.Position(2);jointData.Velocity(1);jointData.Velocity(2)];
    theta_1(i) = jointData.Position(1);
    theta_1_dot(i) = jointData.Velocity(1);
    theta_2(i) = jointData.Position(2);
    theta_2_dot(i) = jointData.Velocity(2);
    
    tau1.Data = -K(1,:)*X;
    tau2.Data = -K(2,:)*X;


    u1(i) = tau1.Data;
    u2(i) = tau2.Data;

    send(j1_effort,tau1);
    send(j2_effort,tau2);
    
    i = i+1;
    
% sample the time, joint state values, and calculated torques here to be
%plotted at the end
end
tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

% disconnect from roscore
%rosshutdown;

% Plot the trajectories

%figure;
subplot(2,3,1);
plot(time,u1)
xlabel('time t in sec');
ylabel('u1');
grid on;

%figure;
subplot(2,3,2);
plot(time,u2)
xlabel('time t in sec');    
ylabel('u2');
grid on;

%figure;
subplot(2,3,3);
plot(time,theta_1);
xlabel('time t in sec');
ylabel('theta_1');
grid on;

%figure;
subplot(2,3,4);
plot(time,theta_1_dot);
xlabel('time t in sec');
ylabel('theta_1_dot');
grid on;

%figure;
subplot(2,3,5);
plot(time,theta_2);
xlabel('time t in sec');
ylabel('theta_2');
grid on;

%figure;
subplot(2,3,6);
plot(time,theta_2_dot);
xlabel('time t in sec');
ylabel('theta_2_dot');
grid on;