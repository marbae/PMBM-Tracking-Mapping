function state = stateTurn(time,L,R,x0,v)
    T_line = abs(L/v(1));
    T_half_circle = pi*abs(R/v(1));
    time = mod(time, 2*T_line+2*T_half_circle);
    if time <= T_line
        x = x0(1) + v(1)*time;
        y = x0(2) + v(2)*time;   
        yaw = x0(3) + v(3)*time;
        vx = v(1);
        vy = v(2);
        yaw_rate = v(3);
    elseif time <= T_line+T_half_circle
        yaw = x0(3) + v(1)/R*(time - T_line);
        yaw_rate = v(3) + v(1)/R;
        x = x0(1) + L + R*cos(yaw - pi/2);
        y = x0(2) + R + R*sin(yaw - pi/2);
        vx = -v(1)*sin(yaw - pi/2);
        vy =  v(1)*cos(yaw - pi/2);            
    elseif time <= 2*T_line+T_half_circle
        x = x0(1) + L - v(1)*(time - T_line - T_half_circle);
        y = x0(2) + 2*R;
        yaw = pi+x0(3);
        vx = -v(1);
        vy = v(2);
        yaw_rate = 0;
    else
        yaw = v(1)/R*(time - 2*T_line - T_half_circle) + pi+x0(3);
        yaw_rate = v(1)/R;
        x = x0(1) + R*cos(yaw - pi/2);
        y = x0(2) + R + R*sin(yaw - pi/2);
        vx = -v(1)*sin(yaw - pi/2);
        vy =  v(1)*cos(yaw - pi/2);
    end
    yaw = ssa(yaw);
    state = [x;y;yaw;vx;vy;yaw_rate];
end