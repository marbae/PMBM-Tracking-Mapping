function f6 = generateSubplot(metrics,GOSPA,type,interval,f6,linestyle)
    Nt = length(metrics);
    numMC = size(metrics,2);
    IOU = zeros(1,Nt);
    Gospa = zeros(1,Nt);
    locErr = zeros(1,Nt);
    missedT = zeros(1,Nt);
    falseT = zeros(1,Nt);
    headingErr = zeros(1,Nt);
    for i = 1:Nt
        Gospa(i) = mean(GOSPA(i,1,:));
        locErr(i) = mean(GOSPA(i,2,:));
        missedT(i) = mean(GOSPA(i,3,:));
        falseT(i) = mean(GOSPA(i,4,:));
        nPresentTargets = numMC;
        for j = 1:numMC
            if(~isempty(metrics{i,j}))
                IOU(i)=IOU(i)+mean(metrics{i,j}.IOU);
                headingErr(i)=headingErr(i)+mean(metrics{i,j}.headingErr);
            else
                nPresentTargets = nPresentTargets-1;
            end
        end
        if(nPresentTargets>0)
            IOU(i) = IOU(i)/nPresentTargets;
            headingErr(i) = headingErr(i)/nPresentTargets;
        end
    end
    figure(f6)
    axis([interval(1) interval(2) 0 inf]);
    hold on
    subplot(3,2,1)
    axis([interval(1) interval(2) 0 inf]);
    hold on
    plot(Gospa,'LineWidth',1,'DisplayName',type,'LineStyle',linestyle)
    title('GOSPA')
    subplot(3,2,2)
    axis([interval(1) interval(2) 0 inf]);
    hold on
    plot(locErr,'LineWidth',1,'DisplayName',type,'LineStyle',linestyle)
    title('Localization Error')
    subplot(3,2,3)
    axis([interval(1) interval(2) 0 inf]);
    hold on
    plot(missedT,'LineWidth',1,'DisplayName',type,'LineStyle',linestyle)
    title('Missed Targets')
    subplot(3,2,4)
    axis([interval(1) interval(2) 0 inf]);
    hold on
    plot(falseT,'LineWidth',1,'DisplayName',type,'LineStyle',linestyle)
    title('False Targets')
    subplot(3,2,5)
    axis([interval(1) interval(2) 0 inf]);
    hold on
    plot(IOU,'LineWidth',1,'DisplayName',type,'LineStyle',linestyle)
    title('IOU')
    subplot(3,2,6)
    axis([interval(1) interval(2) 0 inf]);
    hold on
    plot(headingErr,'LineWidth',1,'DisplayName',type,'LineStyle',linestyle)
    title('Heading Error (rad)')
end

