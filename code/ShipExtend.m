classdef ShipExtend
	properties
        L % length aft (stern) to bow
        W % width port to starboard
        D % distance from bow, where the ship is widest (first)
        S % width of flat stern
        type % type of ship hull
        angles % vector of angles 0<=angles(i)<2*pi
        hull_r % radius values of extend : N_angles x 1
        hull_b % extend coordinates in body frame : 2 x N_angles
    end

    methods
        function obj = ShipExtend(L, W, D, S, type, dangles)
        obj.L = L;
        obj.W = W;
        obj.D = D;
        obj.S = S;
        obj.type = type;
        angles = 0:dangles:2*pi-dangles;
        obj.angles = angles;
        hull_r = getHullRadii(L, W, D, S, type, angles);
        obj.hull_r = hull_r;
        [hull_x, hull_y] = pol2cart(obj.angles(:), hull_r);
        obj.hull_b = [hull_x' hull_x(1);
                      hull_y' hull_y(1)];
        end
	end
end

function hull_r = getHullRadii(L, W, D, S, type, angles)
    % returns radius function values
    % hull_r : N_angles x 1
    Nangles = length(angles);
    if strcmp(type, 'ellipsis')
        x_interpol = linspace(-L/2,L/2,Nangles);
        y_interpol = W/2*sqrt(1 - (x_interpol/L*2).^2);
        [angles_interpol, r_interpol] = cart2pol(x_interpol, y_interpol);
        hull_r = interp1(angles_interpol, r_interpol, abs(ssa(angles)), 'spline');
    elseif strcmp(type, 'rectangle')
        x_interpol = zeros(1,Nangles);
        y_interpol = zeros(1,Nangles);
        ratio = W/L;
        R_x = L/(2*cos(asin(ratio)));
        for i = 1:Nangles
            if((abs(L/2*sin(angles(i)))<W/2) && cos(angles(i))>0)
                x_interpol(i) = L/2;
                y_interpol(i) = L/2*sin(angles(i));
            elseif((abs(L/2*sin(angles(i)))>=W/2) && sin(angles(i))>0)
                x_interpol(i) = R_x*(cos(angles(i)));
                y_interpol(i) = W/2;
            elseif((abs(L/2*sin(angles(i)))<W/2) && cos(angles(i))<0)
                x_interpol(i) = -L/2;
                y_interpol(i) = L/2*sin(angles(i));
            elseif((abs(L/2*sin(angles(i)))>=W/2) && sin(angles(i))<0)
                x_interpol(i) = R_x*(cos(angles(i)));
                y_interpol(i) = -W/2;
            end
        end
        [angles_interpol, r_interpol] = cart2pol(x_interpol, y_interpol);
        hull_r = interp1(angles_interpol, r_interpol, abs(ssa(angles)), 'spline');
    elseif strcmp(type, 'parabola')
        x_interpol = linspace(-L/2,L/2,Nangles);
        y_interpol = W/2 - 2*W/L^2*x_interpol.^2;
        [angles_interpol, r_interpol] = cart2pol(x_interpol, y_interpol);
        hull_r = interp1(angles_interpol, r_interpol, abs(ssa(angles)), 'spline');
    elseif strcmp(type, 'ellipsisFlatStern')
        m = -W^2/4/D^2;
        n = -W^2/2/D;
        p = -(W^2-S^2)/4/(L-D)^2;
        q =  (W^2-S^2)/2/(L-D);
        x1 = linspace(-L/2,(L/2-D),Nangles);
        x1 = x1(1:end-1);
        x2 = linspace((L/2-D),L/2,Nangles);
        x_interpol = [x1, x2];
        y1 = p*(x1+L/2).^2 + q*(x1+L/2) + S^2/4;
        y2 = m*(x2-L/2).^2 + n*(x2-L/2);
        y_interpol = [y1, y2];
        y_interpol = sqrt(y_interpol);
                [angles_interpol, r_interpol] = cart2pol(x_interpol, y_interpol);
        theta_S = atan2(S,-L);
        for n = Nangles:-1:1
            theta = abs(ssa(angles(n)));
            if theta <= theta_S
                hull_r(n) = interp1(angles_interpol, r_interpol, theta, 'spline');
            else
                hull_r(n) = -L/2/cos(theta);
            end
        end
    elseif strcmp(type, 'parabolaFlatStern')
        m = -W/2/D^2;
        n = -W/D;
        p = -(W-S)/2/(L-D)^2;
        q =  (W-S)/(L-D);
        x1 = linspace(-L/2,(L/2-D),Nangles);
        x1 = x1(1:end-1);
        x2 = linspace((L/2-D),L/2,Nangles);
        x_interpol = [x1, x2];
        y1 = p*(x1+L/2).^2 + q*(x1+L/2) + S/2;
        y2 = m*(x2-L/2).^2 + n*(x2-L/2);
        y_interpol = [y1, y2];
        [angles_interpol, r_interpol] = cart2pol(x_interpol, y_interpol);
        theta_S = atan2(S,-L);
        for n = Nangles:-1:1
            theta = abs(ssa(angles(n)));
            if theta <= theta_S
                hull_r(n) = interp1(angles_interpol, r_interpol, theta, 'spline');
            else
                hull_r(n) = -L/2/cos(theta);
            end
        end
    elseif strcmp(type, 'boxEllipticBow')
        x_interpol = linspace(L/2-D,L/2,Nangles);
        y_interpol = W/2*sqrt(1 - ((x_interpol-(L/2-D))/D).^2);
        [angles_interpol, r_interpol] = cart2pol(x_interpol, y_interpol);
        theta1 = atan2(W/2,L/2-D);
        theta2 = atan2(W,-L);
        for n = Nangles:-1:1
            theta = abs(ssa(angles(n)));
            if theta <= theta1
                hull_r(n) = interp1(angles_interpol, r_interpol, theta, 'spline');
            elseif theta <= theta2
                hull_r(n) = W/2/sin(theta);
            else
                hull_r(n) = -L/2/cos(theta);
            end
        end
    elseif strcmp(type, 'boxParabolicBow')
        x_interpol = linspace(L/2-D,L/2,Nangles);
        y_interpol = -W/2/D^2*(x_interpol-(L/2-D)).^2 + W/2;
        [angles_interpol, r_interpol] = cart2pol(x_interpol, y_interpol);
        theta1 = atan2(W/2,L/2-D);
        theta2 = atan2(W,-L);
        for n = Nangles:-1:1
            theta = abs(ssa(angles(n)));
            if theta <= theta1
                hull_r(n) = interp1(angles_interpol, r_interpol, theta, 'spline');
            elseif theta <= theta2
                hull_r(n) = W/2/sin(theta);
            else
                hull_r(n) = -L/2/cos(theta);
            end
        end
    else
        error('Unknown type. Error.')
    end
    hull_r = hull_r(:);
end
