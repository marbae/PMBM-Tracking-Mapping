classdef GP
    properties
        signal_var
        bias_var
        noise_var
        lengthscale
        symmetric %true or false
        theta_test
        Ntest
        theta_plot
        forgetting_factor
        Ktt
        Test2Plot
%         Kposterior
    end
    
    methods
        function obj = GP(hyperparameters, symmetric, n_theta_test, dtheta_plot, forgetting_factor)
            signal_var = hyperparameters(1)^2;
            obj.signal_var = signal_var;
            bias_var = hyperparameters(2)^2;
            obj.bias_var = bias_var;
            noise_var = hyperparameters(3)^2;
            obj.noise_var = noise_var;
            lengthscale = hyperparameters(4)^2;
            obj.lengthscale = lengthscale;
            obj.symmetric = symmetric;
            dtheta_test = 2*pi/n_theta_test;
            if symmetric
            	theta_test = 0:dtheta_test:2*pi-dtheta_test;
            else
                theta_test = 0:dtheta_test:2*pi-dtheta_test;
            end
            obj.theta_test = theta_test;
            obj.Ntest = length(theta_test);
            theta_plot = 0:dtheta_plot:2*pi-dtheta_plot;
            obj.theta_plot = theta_plot;
            obj.forgetting_factor = forgetting_factor;
            obj.Ktt = covarianceMatrix(theta_test, theta_test, ...
                signal_var, bias_var, noise_var, lengthscale,symmetric);
            Kpt = covarianceMatrix(theta_plot, theta_test, ...
                signal_var, bias_var, noise_var, lengthscale,symmetric);
%             Kpp = covarianceMatrix(theta_plot, theta_plot, ...
%                 signal_var, bias_var, noise_var, lengthscale);
            obj.Test2Plot = Kpt/obj.Ktt;
%             obj.Kposterior = Kpp - Kpt*(obj.Ktt\Kpt');
        end
        
        function F = F(obj, dt)
            F = exp(-dt*obj.forgetting_factor)*eye(obj.Ntest);
            %F = eye(obj.Ntest);
        end
        
        function Q = Q(obj, dt)
            Q = (1-exp(-2*dt*obj.forgetting_factor))*obj.Ktt;
%             obj.forgetting_factor = 1-obj.forgetting_factor;
%             Q = (1/obj.forgetting_factor-1)*obj.Ktt;
        end
        
        function K = covarianceMatrix(obj, thetas1, thetas2)
            K = covarianceMatrix(thetas1, thetas2, obj.signal_var, obj.bias_var, obj.noise_var, obj.lengthscale, obj.symmetric);
        end
        
        function dK = dCovarianceMatrix(obj, thetas1, thetas2)
            dK = dCovarianceMatrixDTheta1(thetas1, thetas2, obj.signal_var, obj.bias_var, obj.noise_var, obj.lengthscale, obj.symmetric);
        end
        
        function extendBody = getExtendBody(obj, extendRadiiTest)
            % returns body coordinates of extend
            % the coordinates give a closed curve
            % hulls_r: Ntest x Nt
            % extendBody : 2 x Nplot+1 x Nt
            Nt = size(extendRadiiTest,2);
            for nt=Nt:-1:1
                extendRadiiPlot = obj.Test2Plot*extendRadiiTest(:,nt);
                [extendBodyX, extendBodyY] = pol2cart(obj.theta_plot(:), extendRadiiPlot);
                extendBody(:,:,nt) = [extendBodyX' extendBodyX(1);
                                    extendBodyY' extendBodyY(1)];
            end
        end

        function extendNE = getExtendNE(obj, extendRadiiTest, positions, headings)
            % returns cartesian coordinates of hull in ned frame
            % the coordinates give a closed curve
            % extendRadiiTest: Ntest x Nt
            % positions : 2 x Nt
            % headings: Nt x 1
            % extendNE: 2 x Nplot x Nt
            extendBody = obj.getExtendBody(extendRadiiTest);
            Nt = size(extendBody,3);
            for nt=Nt:-1:1
                extendNE(:,:,nt) = rotate(extendBody(:,:,nt), headings(nt)) + positions(:,nt);
            end
        end

        function [rTestPred, PPred, rTestUpd, PUpd, innovation, S] = regression(obj, rTestEst, PEst, inputRadii, inputAngles, dt)
            F = obj.F(dt);
            Q = obj.Q(dt);
            if isempty(inputRadii)
                rTestUpd = rTestEst;
                PUpd = PEst;
            else
                I = eye(size(F));
                Kii = covarianceMatrix(obj, inputAngles, inputAngles);
                Kit = covarianceMatrix(obj, inputAngles, obj.theta_test);
                Kti = Kit.';
                H = Kit/obj.K_tt;
                R = Kii - Kit*(obj.Ktt\Kti);
                innovation = inputRadii - H*rTestEst;
                S = H*PEst*H.' + R;
                K = PEst*H.'/S;
                rTestUpd = rTestEst + K*innovation;
                PUpd = (I-K*H)*PEst*(I-K*H).' + K*R*K.';
            end
            rTestPred = F*rTestUpd;
            PPred = F*PUpd*F.' + Q;
        end
    end
end

function K = covarianceMatrix(thetas1, thetas2, s_var, b_var, n_var, l, symmetric)
   N1 = length(thetas1);
   N2 = length(thetas2);
   K = zeros(N1,N2);
   for n1 = 1:N1
       for n2 = 1:N2
           K(n1,n2) = covarianceFunction(thetas1(n1), thetas2(n2), s_var, b_var, n_var, l, symmetric);
       end
   end
end

function dKdTheta1 = dCovarianceMatrixDTheta1(thetas1, thetas2, s_var, b_var, n_var, l , symmetric)
   N1 = length(thetas1);
   N2 = length(thetas2);
   dKdTheta1 = zeros(N1,N2);
   for n1 = 1:N1
       for n2 = 1:N2
           dKdTheta1(n1,n2) = dCovarianceFunctionDTheta1(thetas1(n1), thetas2(n2), s_var, b_var, n_var, l, symmetric);
       end
   end
end

function k = covarianceFunction(theta1, theta2, s_var, b_var, n_var, l, symmetric)            
     if symmetric
         %k = s_var*exp(-(1/(2*l^2))*(sin(theta1-theta2))^2);
         k = s_var*exp(-(1/(2*l^2))*(abs(ssa(theta1))-abs(ssa(theta2)))^2);
%         k = s_var*exp(-(2/l^2)*sin(abs(theta1-theta2)/2)^2)...
%             +s_var*exp(-(2/l^2)*sin(abs(theta1+theta2)/2)^2);
     else
         k = s_var*exp(-(2/l^2)*sin(((theta1-theta2))/2)^2);
     end
    k = k + b_var;
    if abs(theta1-theta2) < 1e-16
        k = k + n_var;
    end
end

function dkdtheta1 = dCovarianceFunctionDTheta1(theta1, theta2, s_var, b_var, n_var, l,symmetric)            
    %k = covarianceFunction(theta1, theta2, s_var, b_var, n_var, l);
    if(symmetric)
        k = s_var*exp(-(1/(2*l^2))*(abs(ssa(theta1))-abs(ssa(theta2)))^2);
        %k = s_var*exp(-(1/(2*l^2))*(sin(theta1-theta2))^2);
        dkdtheta1 = k*(-1/l^2)*(abs(ssa(theta1))-abs(ssa(theta2)))*sign(ssa(theta1));
        %dkdtheta1 = k*(-1/(2*l^2))*(sin(2*(theta1-theta2)));
    else
        k = s_var*exp(-(2/l^2)*sin(((theta1-theta2))/2)^2);    
        dkdtheta1 = k*(-1/(l^2))*(sin(theta1-theta2));
    end

end
