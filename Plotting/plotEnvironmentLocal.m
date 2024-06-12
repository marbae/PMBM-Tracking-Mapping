function [] = plotEnvironmentLocal(egoState,environment)
    egoHeading = deg2rad(egoState.heading);
    egoN = egoState.N;
    egoE = egoState.E;

    for i=1:length(environment)
        feature = environment{i};
        feature(:,1) = feature(:,1)-egoN;
        feature(:,2) = feature(:,2)-egoE;
        r = min(sqrt(feature(:,1).^2+feature(:,2).^2));
        if(r<100)
            rotFeature = feature;
            rotFeature(:,1) = feature(:,1)*cos(egoHeading)-feature(:,2)*sin(egoHeading);
            rotFeature(:,2) = feature(:,1)*sin(egoHeading)+feature(:,2)*cos(egoHeading);
            plot(rotFeature(:,2),rotFeature(:,1))
        end
    end
end

    