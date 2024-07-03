clc;clear all; close all
rng(123)
N = 2000;
t_squ=1:N;
samp_freq_US = 25; % Sampling Frequency for Ultrasound (Hz) 

%% ===== Simulating the moving trajectory of needle tip ===========
% set input velocity of needle
v_bias = 0.001*sin(t_squ/250*pi); % Velocity fluctuation term
v_input = 0.08*[(1:100)/100, ones(1,400), 1.5*ones(1,1000) ones(1,500)]+ v_bias;

% set x y z velocity of needle
vy_needle = v_input.*(1600*[t_squ(1:1200).^2/N^3 t_squ(1200)^2/N^3*ones(1,N-1200)]);
vz_needle = v_input.*(700*[t_squ(1:1500)/N^2 t_squ(1500)/N^2*ones(1,N-1500)]);
vx_needle = sqrt(v_input.^2 - vy_needle.^2 - vz_needle.^2);

% get x y z trajectories (true) of needle
x_needle=cumsum(vx_needle);
y_needle=cumsum(vy_needle);
z_needle=cumsum(vz_needle);

% % ----- show the planned trajectories ----
% figure('Position',[100,100,1200,800])
% subplot(3,2,1)
% plot(x_needle);
% title('Moving curve of the needle tip in the x-direction');
% xlabel('Samples (1/25 s)');
% ylabel('x (mm)');
% subplot(3,2,3)
% plot(y_needle);
% title('Moving curve of the needle tip in the y-direction');
% xlabel('Samples (1/25 s)');
% ylabel('y (mm)');
% subplot(3,2,5)
% plot(z_needle);
% title('Moving curve of the needle tip in the z-direction');
% xlabel('Samples (1/25 s)');
% ylabel('z (mm)');
% 
% subplot(3,2,2)
% plot(v_input*samp_freq_US,'r'); hold on
% plot(vx_needle*samp_freq_US,'g');
% plot(vy_needle*samp_freq_US,'b');
% plot(vz_needle*samp_freq_US,'y');
% legend('v_{in}','v_{x}', 'v_{y}', 'v_{z}');
% title('Input velocity (v_{in}) and tip velocities in x, y, z directions');
% xlabel('Samples (1/25 s)');
% ylabel('velocity (mm/s)');
% ylim([-0.2, 3.3]);
% 
% subplot(3,2,[4 6])
% plot3(x_needle,y_needle, z_needle);hold on;
% xlim([0, 200]);
% ylim([-100, 100]);
% zlim([-100, 100]);
% title('Moving trajectory of the needle tip in 3D space');
% xlabel('x (mm)');
% ylabel('y (mm)');
% zlabel('z (mm)');
% grid on

%% == Initialization for tracking ==
p_measure = zeros(3,N); % Postion measurments obtained by Ultrasound
state_squ = zeros(5,N); % State vectors for Kalman filtering

% state_squ(:,1) = [0;0;0;pi/2;pi/2]; % for EKF_PositionAngles
%state_squ(:,1) = [0;0;0;0;0]; % for EKF_PositionVelocityRatio
state_squ(:,1) = [0;0;0;vy_needle(1);vz_needle(1)]; % for EKF_PositionVelocity.

Pe = 0; % initialize the estimate error covariance
vx_EKF = zeros(1,N); % x velocity estimated by EKF
vy_EKF = zeros(1,N); % y velocity estimated by EKF

y_US = zeros(1,N); % y position of the ultrasound probe
x_US = zeros(1,N); % x position of the ultrasound probe
x_US(1:2) = x_needle(1);
vx_US = zeros(1,N); % x velocity of the ultrasound probe

v_e_sum = 0; % velocity error accumulation for PI controlling of x velocity
xIds = 2;
behandId = 1;
aheadId = 0;

%% == Run tracker ==
for i=2:2000    
    % --- get noisy measurements based on the x position of ultrasound.(x_US)
    if x_US(i) > x_needle(i) % Ultrasound probe running ahead of the needle
        R=0;
        p_measure(:,i) = p_measure(:,i-1); % No measurement, so use the previous measurement.
        aheadId = i;
        Kvx = 0.80^(aheadId-behandId); % Our strategy
%         Kvx = 0.5; % 2013 ICRA paper
    else % Ultrasound probe behind the needle tip
        R=[0.5 0 0;
           0 100 0; 
           0 0 100];
        xIds = (xIds-2) + find(x_needle((xIds-1):i)>x_US(i),1);
        ratio_x = (x_US(i)-x_needle(xIds-1))/(x_needle(xIds)-x_needle(xIds-1));
        y_measure = y_needle(xIds-1) + (y_needle(xIds)-y_needle(xIds-1))*ratio_x + 4*rand(1)-2;
        z_measure = z_needle(xIds-1) + (z_needle(xIds)-z_needle(xIds-1))*ratio_x + 4*rand(1)-2;      
        p_measure(:,i) = [x_US(i); y_measure; z_measure]; % Simulation of noisy measurements
        behandId = i;
        Kvx = 1.01^(behandId-aheadId); % Our strategy
%         Kvx = 1.05; % 2013 ICRA paper
    end
    % --- get noisy measurements based on the current X position of ultrasound.(x_US)
    
    % --- EKF based on p_measure ---   
%     [state_squ(:,i), Pe, vx_EKF(i), vy_EKF(i)] = EKF_PositionAngles(state_squ(:,i-1), Pe, v_input(i), p_measure(:,i), R);
%     [state_squ(:,i), Pe, vx_EKF(i), vy_EKF(i)] = EKF_PositionVelocityRatio(state_squ(:,i-1), Pe, v_input(i), p_measure(:,i), R);
    [state_squ(:,i), Pe, vx_EKF(i), vy_EKF(i)] = EKF_PositionVelocity(state_squ(:,i-1), Pe, v_input(i), p_measure(:,i), R);
    
    %====== Execute PD poisiton control in the Y direction (lasts for one sample interval)
    [y_US_sub, ~] = PDPositionControl(y_US(i), state_squ(2,i));
    y_US(i+1) = y_US_sub(end); % results of PD control,i.e. the ultrasound probe moves to Y_US
    
    %====== Execute PI velcotiy control in the X direction 
    [x_v_US_sub, v_e_sum, delta_s] = PIVelocityControl(vx_US(i), v_e_sum, Kvx*vx_EKF(i));
    vx_US(i+1) = x_v_US_sub(end); % x velocity of ultrasound after PD control
    x_US(i+1) = x_US(i)+delta_s; % x position of ultrasound probe after PD control
end

%% === Evaluating the results of tracking ==
pos_err = state_squ(1:3,:)-[x_needle;y_needle;z_needle];
filter_rmse = sqrt(sum(sum((pos_err(:,3:N)).^2))/(N-2));  % RMSE of tip position estimation with Kalman filter
filter_vx_rmse = sqrt(sum((vx_EKF(3:N)-vx_needle(3:N)).^2)/(N-2));

x_US_rmse = sqrt(sum((x_US(3:N)-x_needle(3:N)).^2)/(N-2)); % x-direction RMSE between ultrasound and needle tip
y_US_rmse = sqrt(sum((y_US(3:N)-y_needle(3:N)).^2)/(N-2)); % y-direction RMSE between ultrasound and needle tip

isAhead = x_US(1:N)>x_needle(1,:);
lead_lag_ratio = sum(isAhead)/(N-sum(isAhead));

%----- show results ------
figure('Position',[100,100,1200,800])
subplot(2,2,1)
plot(p_measure(1,:),'y'); hold on;
plot(x_needle,'r'); 
plot(state_squ(1,:),'g'); 
plot(x_US(1:N),'b'); 
title('Position curves in the x-direction');
h1=legend('$z_{x}$','$p_{x}$', '$ \hat p_{x}$','$p_{ux}$','latex','Location','southeast' )
xlabel('Samples (1/25 s)');
ylabel('x (mm)');
set(h1,'Interpreter','latex')

subplot(2,2,3)
plot(p_measure(2,:),'y'); hold on;
plot(y_needle,'r');
plot(state_squ(2,:),'g'); 
plot(y_US,'b'); 
title('Position curves in the y-direction');
h2=legend('$z_{y}$', '$p_{y}$', '$ \hat p_{y}$','$p_{uy}$','latex','Location','southeast' )
xlabel('Samples (1/25 s)');
ylabel('y (mm)');
set(h2,'Interpreter','latex')

% 绘制y方向的位置曲线
subplot(2,2,3)
plot(1:2000, p_measure(2,1:2000),'y'); hold on;
plot(1:2000, y_needle(1:2000),'r');
plot(1:2000, state_squ(2,1:2000),'g'); 
plot(1:2000, y_US(1:2000),'b'); 
title('Position curves in the y-direction');
h3=legend('$z_{y}$', '$p_{y}$', '$ \hat p_{y}$','$p_{uy}$','Location','southeast' )
xlabel('Samples (1/25 s)');
ylabel('y (mm)');
xlim([1 2000]); % 设置横坐标范围
set(h3,'Interpreter','latex')


% 绘制速度曲线
subplot(2,2,2)
plot(1:2000, vx_needle*samp_freq_US); hold on
plot(1:2000, vx_EKF*samp_freq_US);
plot(1:2000, vy_needle*samp_freq_US);
plot(1:2000, vy_EKF*samp_freq_US);
ylim([-1.2, 3.25]);
title('Velocity curves in the x- and y-directions');
h4=legend('$v_{x}$', '$ \hat v_{x}$', '$v_{y}$','$ \hat v_{y}$','latex','Location','southeast' )
xlabel('Samples (1/25 s)');
ylabel('Velocity (mm/s)');
xlim([1 2000]); % 设置横坐标范围
set(h4,'Interpreter','latex')

subplot(2,2,4)
plot3(p_measure(1,:),p_measure(2,:),p_measure(3,:),'y');hold on
plot3(x_needle,y_needle,z_needle,'r');
plot3(state_squ(1,:),state_squ(2,:),state_squ(3,:),'g'); 
xlim([0, 200]);
ylim([-100, 100]);
zlim([-100, 100]);
grid on
title('3D trajectories');
legend('measures', 'needle', 'EKF','Location','northeast' )
xlabel('x/mm')
ylabel('y/mm')
zlabel('z/mm')
