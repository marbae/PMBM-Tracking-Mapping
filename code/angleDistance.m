function distance = angleDistance(angle1,angle2)
    angle1 = ssa(angle1);
    angle2 = ssa(angle2);
    distance = abs(angle1-angle2);
    if(distance>pi)
        distance = 2*pi-distance;
    end
end