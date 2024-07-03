function [state, P, vx_filted, vy_filted] = EKF_PositionVelocityRatio(state_pre, P_pre, v_input, p_measure, R)

    % Calculate a priori estimate of current state based on previous state
    % state_pre - previous state
    cos_alpha_pre = sqrt(1-state_pre(4)^2-state_pre(5)^2);
    state_ = state_pre + [v_input*cos_alpha_pre; 
                          v_input*state_pre(4); 
                          v_input*state_pre(5); 
                          0;
                          0];  % Estimating a priori state using uniform velocity model 
                      
    % Calculate the linearized state transfer matrix (Jacobi) at the current moment
    A=[1 0 0 -v_input*state_pre(4)/cos_alpha_pre  -v_input*state_pre(5)/cos_alpha_pre;
       0 1 0 v_input 0;
       0 0 1  0 v_input;
       0 0 0 1 0;
       0 0 0 0 1];
   
    % set the observation matrix
    H=[1 0 0 0 0;
        0 1 0 0 0;
        0 0 1 0 0]; 
    % set the process noise covariance matrix
    Q=[5e-6 0 0 0 0;
        0 10e-6 0 0 0;
        0 0 10e-6 0 0;
        0 0 0 5e-4 0;
        0 0 0 0 5e-4];
   
    % Kalman filtering
    [state, P] = Kalman_filter3(0, P_pre, p_measure, A, H, Q, R, state_);
    if (state(4)^2+state(5)^2)>1
        vx_filted = 0; 
    else
        vx_filted = v_input*sqrt(1-state(4)^2-state(5)^2);
    end
    vy_filted = v_input*state(4);