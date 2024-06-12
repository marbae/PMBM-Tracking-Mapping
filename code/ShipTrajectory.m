classdef ShipTrajectory
    properties
        times %vector of time points: 1 x Nt
%         x0 % initial state: 6 x 1, [N;E;psi;v_N;v_E;r] or 4x1 [N;E;v_N;v_E]        
%         q % q for cv2 or cv3 model        
        extendAngles % Nangles x 1
        extendRadii % Nangles x 1
        extendBody % 2 x Nangles+1
        xKin % kinematic states: 6 x Nt
        extendNE % ship/hull/extend trajectory: 2 x Nangles+1 x Nt
        lowerextendNE % ship/hull/extend trajectory: 2 x Nangles+1 x Nt
        upperextendNE % ship/hull/extend trajectory: 2 x Nangles+1 x Nt
        motionType % String;
    end
    
    methods
        function obj = ShipTrajectory(times, q, x0, extend, motionType)
            obj.times = times;
            obj.motionType = motionType;
%             obj.x0 = x0;
%             obj.q = q;
            obj.extendAngles = extend.angles';
            obj.extendRadii = extend.hull_r;
            extendBody = extend.hull_b;
            obj.extendBody = extendBody;
            xKin = getKinematicStates(times, x0, q, motionType);
            obj.xKin = xKin;
            obj.extendNE = getExtendNE(extendBody, xKin(1:2,:), xKin(3,:)');
            obj.lowerextendNE = extend.hull_b;
            obj.upperextendNE = extend.hull_b;
        end
    end
end

function xKin = getKinematicStates(times, x0, q,motionType)
    % xKin: 6 x Nt
    Nt = length(times);
    dT = diff(times);            
    xKin = zeros(6,Nt);
    cv = CV(q);
    xt = x0;
    L=40;
    R=20;
    % save initial state
    if cv.NDoF == 2
        xKin([1 2 4 5], 1) = xt;
    elseif cv.NDoF == 3
        xKin(1:6, 1) = xt;
    end
    % simulate state and save
    for nt = 2:Nt
        dt = dT(nt-1);
        if(strcmp(motionType,'random'))
            xt = cv.F(dt)*xt + ...
            chol(cv.Q(dt))'*randn(2*cv.NDoF,1);
        end
        if(strcmp(motionType,'constant'))
            xt = cv.F(dt)*xt;
        end
        if(strcmp(motionType,'constantTurn'))
            xt=stateTurn(nt*dt,L,R,x0(1:3),x0(4:6));
        end
        if cv.NDoF == 2
            xKin([1 2 4 5], nt) = xt;
        elseif cv.NDoF == 3
            xKin(1:6, nt) = xt;
        end
    end
    % calculate psi and r when cv2 used
    if cv.NDoF == 2
        xKin(3,:) = atan2(xKin(5,:),xKin(4,:))';
        if Nt == 1
            xKin(6,:) = NaN;
        elseif Nt == 2
            xKin(6,:) = (xKin(3,2)-xKin(3,1))/dT(1);
        else
            for nt = 1:Nt
                if nt == 1
                    xKin(6,nt) = (xKin(3,2)-xKin(3,1))/dT(1);
                elseif nt == Nt
                    xKin(6,nt) = (xKin(3,Nt)-xKin(3,Nt-1))/dT(end);
                else
                    xKin(6,nt) = (xKin(3,nt+1)-xKin(3,nt-1))/(dT(nt)+dT(nt-1));
                end
            end
        end
    end
end

function extendNE = getExtendNE(extendBody, positions, headings)
    % returns cartesian coordinates of extend trajectory in ned frame
    % extendBody: 2 x Nangles+1
    % positions : 2 x Nt
    % headings: Nt x 1
    % hulls_n: 2 x Nangles+1 x Nt
    Nt = length(headings);
    for nt=Nt:-1:1
     extendNE(:,:,nt) = rotate(extendBody, headings(nt)) + positions(:,nt);
    end
end
