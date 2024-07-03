function [s_real, v_real] = PDPositionControl(s0, s_obj, Kp, Kd, mess, b, K)
    if nargin<=3
        Kp = 30;
        Kd = 50;
        mess = 150;
        b = 20;
        K = 300; %Number of control cycles
    end

    v_real=zeros(1,K);
    s_real=zeros(1,K);
    s_sensor = s0 + 0.0000008*(rand(1)-0.5);
    v_sensor = 0;
    for k=1:K
        v_real(k)=(Kp*(s_obj-s_sensor)+(mess-Kd)*v_sensor)/(b+mess);
        if k==1
            s_real(k)=s0+v_real(k);
        else
            s_real(k)=s_real(k-1)+v_real(k);
        end
        
        v_sensor = s_real(k) + 0.0000008*(rand(1)-0.5) - s_sensor;
        s_sensor = s_sensor+v_sensor;
    end
