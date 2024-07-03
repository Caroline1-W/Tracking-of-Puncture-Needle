function [state, P, vx_filted, vy_filted] = EKF_PositionVelocity(state_pre, P_pre, v_input, p_measure, R)
%   Calculate a priori estimate of current state based on previous state
%     state_pre = state_squ(:,i-1); % previous state
    vx_pre = sqrt(v_input^2-state_pre(4)^2-state_pre(5)^2);
    state_ = state_pre + [vx_pre; 
                          state_pre(4); 
                          state_pre(5); 
                          0;
                          0];  % Estimating a priori state using uniform velocity model 
                      
    % Calculate the linearized state transfer matrix (Jacobi) at the current moment
    A=[1 0 0 -state_pre(4)/vx_pre  -state_pre(5)/vx_pre;
       0 1 0 1 0;
       0 0 1 0 1;
       0 0 0 1 0;
       0 0 0 0 1];
   % set the observation matrix
    H=[1 0 0 0 0;
       0 1 0 0 0;
       0 0 1 0 0]; 
    % set the process noise covariance matrix
    Q=[10e-6 0 0 0 0;
       0 10e-6 0 0 0;
       0 0 10e-6 0 0;
       0 0 0 5e-6 0;
       0 0 0 0 5e-6]; % Q for EKF_Position_Velocity

    % Kalman filtering
    [state, P] = Kalman_filter3(0, P_pre, p_measure, A, H, Q, R, state_);
    if (state(4)^2+state(5)^2)>v_input^2
        vx_filted = 0; 
    else
        vx_filted = sqrt(v_input^2-state(4)^2-state(5)^2);
    end
    vy_filted = state(4);
end