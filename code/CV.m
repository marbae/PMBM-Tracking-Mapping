classdef CV
    % Constant Velocity model of N DoF
    % States x = [x1;...;xN;v1;...;vN]
    % Noise strength q = [q1;...;qN] unit m^2/s^2
    % NDoF Number Degrees of Freedom
    properties
        q
        NDoF
    end

    methods
        function obj = CV(q)
           obj.q = q;
           obj.NDoF = length(q);
        end

        function F = F(obj, dt)
            F = eye(2*obj.NDoF) + diag(dt * ones(obj.NDoF,1), obj.NDoF);
        end

        function Q = Q(obj, dt)
            Q = kron([dt^2/3 dt/2; dt/2 1],diag(obj.q));
        end
    end
end
