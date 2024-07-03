clc;clear all; close all
M = 1;
filter_rmse_ave = 0;
filter_vx_rmse_ave = 0;
x_US_rmse_ave = 0;
y_US_rmse_ave = 0;
lead_lag_ratio_ave = 0;

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

% ----- show the planned trajectories ----
figure('Position',[100,100,1200,800])
subplot(3,2,1)
plot(x_needle);
title('Moving curve of the needle tip in the x-direction');
xlabel('Samples (1/25 s)');
ylabel('x (mm)');
subplot(3,2,3)
plot(y_needle);
title('Moving curve of the needle tip in the y-direction');
xlabel('Samples (1/25 s)');
ylabel('y (mm)');
subplot(3,2,5)
plot(z_needle);
title('Moving curve of the needle tip in the z-direction');
xlabel('Samples (1/25 s)');
ylabel('z (mm)');

subplot(3,2,2)
plot(v_input*samp_freq_US,'r'); hold on
plot(vx_needle*samp_freq_US,'g');
plot(vy_needle*samp_freq_US,'b');
plot(vz_needle*samp_freq_US,'y');
legend('v_{in}','v_{x}', 'v_{y}', 'v_{z}');
title('Input velocity (v_{in}) and tip velocities in x, y, z directions');
xlabel('Samples (1/25 s)');
ylabel('velocity (mm/s)');
ylim([-0.2, 3.3]);

subplot(3,2,[4 6])
plot3(x_needle,y_needle, z_needle);hold on;
xlim([0, 200]);
ylim([-100, 100]);
zlim([-100, 100]);
title('Moving trajectory of the needle tip in 3D space');
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
grid on

%% Repeat M times for validations
for k = 1:M
    %=======
    tracking
    %=======
    
    disp(['No ', num2str(k), ': ', num2str([filter_rmse, x_US_rmse, y_US_rmse, lead_lag_ratio])]);
    filter_rmse_ave = filter_rmse_ave + filter_rmse;
    filter_vx_rmse_ave = filter_vx_rmse_ave + filter_vx_rmse;
    x_US_rmse_ave = x_US_rmse_ave + x_US_rmse;
    y_US_rmse_ave = y_US_rmse_ave + y_US_rmse;
    lead_lag_ratio_ave = lead_lag_ratio_ave + lead_lag_ratio;
end
filter_rmse_ave/M
filter_vx_rmse_ave*samp_freq_US/M
x_US_rmse_ave/M
y_US_rmse_ave/M
lead_lag_ratio_ave/M