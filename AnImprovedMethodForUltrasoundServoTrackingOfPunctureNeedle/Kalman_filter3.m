function [x, P] = Kalman_filter3(x_pre, P_pre, z, A, H, Q, R, x_, W, V)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   x(k) = A*x(k-1)+W*w(k-1) % the process equation     %
%   z(k) = H*x(k)+V*v(k) % the measurment equation      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if input x == 0 denotes Extended KF

% A: the system transfer matrix (linearization for EKF)
% Q: the process noise covariance matrix
% H: the observation matrix (linearization for EKF)
% R: the measurement noise covriance matrix（先验或后验状�?�估计）
% x_ or x: a priori or posteriori state estimate（先验或后验估计误差协方差）
% P_ or P: a priori or posteriori estimate error covariance（先验或后验估计误差协方差）
% K: gain (blending factor) 
% r: innovation or residual r=(z-z_) where z_ or z: predicted or actual measurement
% W,V: Jacobian matrices for Extended Kalman Filter

if nargin<=8
    W = 1;
    V = 1;
end

% Predict steps:
if x_pre~=0
    x_ = A*x_pre;
end
P_ = A*P_pre*A'+W*Q*W'; 

% if No measurment
if R==0                 
    x = x_;
    P = P_;
    return;
end

% Update steps:
r=z-H*x_ ;%计算新息
% (H*P_*H'+V*R*V')
K = P_*H'/(H*P_*H'+V*R*V');
% x_
x = x_+K*r;
P = P_-K*H*P_;