function []  = plotTimeLapseLocal(egoState,ShipTrajectoriesTrue,environment,ZTrue,mapVis,model,doVideo,doLive,plotRate,movieName)
    % Plots trajectory of 1 ship
    % ShipTrajectoryTrue
    % ZTrue
    % ShipTrajectoryEst
    % NELimits
    % options
    % colors
    NELimits = [-100 100 -100 100];
    colorConfidence = [245,222,179]/255; % wheat
    colorShipTrue =   [0, 0.4470, 0.7410; % blue
                   0.8500, 0.3250, 0.0980; % red
                   0.3010, 0.7450, 0.9330]; % green
    NcolorShipTrue = size(colorShipTrue,1);
    colorShipEst =   [0.9290, 0.6940, 0.1250; % orange
                   0.4940, 0.1840, 0.5560; % yellow
                   0.4660, 0.6740, 0.1880; % lilla
                   0.6350, 0.0780, 0.1840]; % ochre
    NcolorShipEst = size(colorShipEst,1);
    sizeCenterShip = 25;
    sizeLidar = 25;
    colorLidar = 'k';
    sizeMeasurement = 5;
    colorTrueMeasurement = 'k';
    % variables
    NShip = length(ShipTrajectoriesTrue);
    if isempty(environment)
        plotEnv = false;
    else
        plotEnv = true;
    end
    if isempty(ZTrue)
        plotZ = false;
    else
        plotZ = true;
    end
    Nt = length(egoState); %TODO
%     if isempty(ShipTrajectoriesEst)
%         plotEst = false;
%         t_birth=1;
%         Nt = length(ZTrue);
%     else
%         plotEst = true;
%         finalTrajectories = ShipTrajectoriesEst{Nt};
%         NShipEst = length(ShipTrajectoriesEst{Nt});
%         t_birth = finalTrajectories(1).t_birth;
%     end
    if isempty(ShipTrajectoriesTrue)
        plotGT = false;
    else
        plotGT = true;
        if(isa(ShipTrajectoriesTrue,'ShipTrajectoryTest'))
            ShipTrajectoriesTruecell = cell(1,1);
            ShipTrajectoriesTruecell{1} = ShipTrajectoriesTrue;
            ShipTrajectoriesTrue=ShipTrajectoriesTruecell;
        end
    end
        
    % timelapse
    if doVideo
        v = VideoWriter(movieName,'MPEG-4');
        %v.FileFormat = 'mp4';
        v.FrameRate = 10;
        open(v)
    end
    figure
    fig = gcf;
    fig.WindowState = 'maximized';
%     tic 
    for nt=1:Nt
        tic
        clf(fig)
        hold on
%         if plotEst
%             for ns=1:NShipEst
%                 t_birth = finalTrajectories(ns).t_birth; 
%                 if(finalTrajectories(ns).t_birth <= nt && finalTrajectories(ns).t_death >= nt)
%                 plotEstimatedShip(finalTrajectories(ns).xKin{nt-t_birth+1},...
%                     finalTrajectories(ns).extendNE{nt-t_birth+1},...
%                     finalTrajectories(ns).upperextendNE{nt-t_birth+1},...
%                     finalTrajectories(ns).lowerextendNE{nt-t_birth+1},...
%                     sizeCenterShip,...
%                     colorShipEst(mod(ns,NcolorShipEst)+1,:))
%                 end
%             end
%         end
        if plotGT
            for ns=1:NShip
                if(isa(ShipTrajectoriesTrue{ns},'ShipTrajectoryTest'))
                    plotShipLocal(ShipTrajectoriesTrue{ns},...
                        sizeCenterShip,...
                        colorShipTrue(mod(ns,NcolorShipTrue)+1,:),egoState{nt},nt)
                else
                    scatter(ShipTrajectoriesTrue{ns}(2,nt),ShipTrajectoriesTrue{ns}(1,nt),sizeCenterShip,colorShipTrue(mod(ns,NcolorShipTrue)+1,:),'filled')
                end
            end
        end
        if plotEnv
            plotEnvironmentLocal(egoState{nt},environment)
        end
         plotMapLocal(egoState{nt},mapVis.allOccMap(:,:,nt),mapVis.allN(:,:,nt),mapVis.allE(:,:,nt),model)
        if plotZ
           plotLidarSensor(sizeLidar,colorLidar)
           plotMeasurements(ZTrue{nt},sizeMeasurement,colorTrueMeasurement)
        end
        
        axis equal
        axis(NELimits)
        xlabel('E [m]')
        ylabel('N [m]')
        title(movieName, "k= "+string(nt))
        hold off
        if doLive
            drawnow
            if nt == 1
                display('Hit a key to start timelapse')
                pause(plotRate)
                tic
            end
            pause(plotRate - toc)
        end
        if doVideo
            frame = getframe(gcf);
            writeVideo(v,frame)
        end
    end
    if doVideo
      close(v)
    end
end

