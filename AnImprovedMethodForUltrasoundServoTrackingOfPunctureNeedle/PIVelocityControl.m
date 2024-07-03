function [v_real, e_sum, delta_s] = PIVelocityControl(v0, e_sum, v_obj, Kp, Ki, mess, b, K)
    if nargin<=3
        Kp = 100;
        Ki = 20;
        mess = 150;
        b = 20;
        K = 300; %Number of control cycles
    end
    delta_s = 0;
    v_obj = v_obj/K;
    v_real=zeros(1,K);
    v_measure = v0/K + 0.0000008*(rand(1)-0.5);
    for k=1:K
        e_k = v_obj-v_measure;
        e_sum = e_sum + e_k;
        v_real(k)=(Kp*e_k + Ki*e_sum + mess*v_measure)/(b+mess);
        v_measure = v_real(k) + 0.0000016*(rand(1)-0.5);
        delta_s = delta_s+v_real(k);
    end
    v_real = v_real*K;
