close all
axis([-70 100 -100 100])
hold on
stop=600;
plot(shipTrajectory1.xKin(2,:),shipTrajectory1.xKin(1,:),'LineWidth',1)
plot(shipTrajectory2.xKin(2,1:5:stop),shipTrajectory2.xKin(1,1:5:stop),'.','LineWidth',1,MarkerSize=5)
plot(shipTrajectory1.extendNE(2,:,20),shipTrajectory1.extendNE(1,:,20),'LineWidth',2)
plot(shipTrajectory2.extendNE(2,:,1),shipTrajectory2.extendNE(1,:,1),'LineWidth',2)
scatter(0,0,25,'k','filled')
grid on
%% 