function globalPos = convertLocalToGlobal(ego,localPos)
    theta = -deg2rad(ego.heading);
    globalPos(1,:) = localPos(1,:)*cos(theta)-localPos(2,:)*sin(theta);
    globalPos(2,:) = localPos(1,:)*sin(theta)+localPos(2,:)*cos(theta);
    globalPos(1,:) = globalPos(1,:)+ego.N;
    globalPos(2,:) = globalPos(2,:)+ego.E;
end