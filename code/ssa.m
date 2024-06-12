function theta_1 = ssa(theta_0)
    % Smallest Signed Angle
    % -pi<theta_1<=pi
    theta_1 = pi - mod(pi - theta_0, 2*pi);
end