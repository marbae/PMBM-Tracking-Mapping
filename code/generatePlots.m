function [] = generatePlots(metrics,model,targets,numMC)
    close all
    Nt = length(metrics);
    f1 = figure('Name','Intersect over Union');
    f2 = figure('Name','RMSE');
    f3 = figure('Name','Heading error');
    f4 = figure('Name','RMSEvel');
    f5 = figure('Name','NEES');
    f6 = figure('Name','NEESposvel');
    f7 = figure('Name','NEESheading');
    f8 = figure('Name','Angular Velocity');
    for iTar = 1:targets
        xAx = [];
        IOU = [];
        RMSE = [];
        headingErr = [];
        RMSEvel = [];
        NEES = [];
        NEESposvel = [];
        NEESheading = [];
        angvel = [];
        for i = 1:Nt
            if(~isempty(metrics{i,numMC}))
                whichTarget = metrics{i,numMC}.iGT==iTar;
                if(any(whichTarget))
                    xAx=[xAx i];
                    IOU = [IOU metrics{i,numMC}.IOU(whichTarget)];
                    RMSE = [RMSE metrics{i,numMC}.RMSE(whichTarget)];
                    headingErr = [headingErr metrics{i,numMC}.headingErr(whichTarget)];
                    RMSEvel = [RMSEvel metrics{i,numMC}.RMSEvel(whichTarget)];
                    NEES = [NEES metrics{i,numMC}.NEES(whichTarget)];
                    NEESposvel = [NEESposvel metrics{i,numMC}.NEESposvel(whichTarget)];
                    NEESheading = [NEESheading metrics{i,numMC}.NEESheading(whichTarget)];
                    angvel = [angvel metrics{i,numMC}.angvel(whichTarget)];
                    
                end
            end
        end
        figure(f1)
        plot(xAx,IOU);
        hold on
        figure(f2)
        plot(xAx(1:end-1),RMSE(1:end-1));
        hold on
        figure(f3)
        plot(xAx,headingErr);
        hold on
        figure(f4)
        plot(xAx(1:end-1),RMSEvel(1:end-1));
        hold on
        figure(f5)
        plot(xAx(1:end-1),log(NEES(1:end-1)));
        hold on
        figure(f6)
        plot(xAx(1:end-1),log(NEESposvel(1:end-1)));
        hold on
        figure(f7)
        plot(xAx(1:end-1),log(NEESheading(1:end-1)));
        hold on
        figure(f8)
        plot(xAx(1:end-1),angvel(1:end-1))
        hold on

    end
    figure(f5)
    upper=chi2inv(0.975,model.Nx)*ones(1,Nt+1);
    lower=chi2inv(0.025,model.Nx)*ones(1,Nt+1);
    plot(log(upper))
    hold on
    plot(log(lower))
    figure(f6)
    upper=chi2inv(0.975,4)*ones(1,Nt+1);
    lower=chi2inv(0.025,4)*ones(1,Nt+1);
    plot(log(upper))
    hold on
    plot(log(lower))
    figure(f7)
    upper=chi2inv(0.975,2)*ones(1,Nt+1);
    lower=chi2inv(0.025,2)*ones(1,Nt+1);
    plot(log(upper))
    hold on
    plot(log(lower))
end