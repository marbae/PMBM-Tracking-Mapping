function [mean,lower,upper] = getExtendCoordinates(GPState,model)
%GETEXTENDCOORDINATES Gets NE coordinates from a curve and a shape (GPState
%   Detailed explanation goes here
    xUpd = GPState.state;
    P = GPState.cov;
    extendStart = model.Nkin+1;
    stdState = sqrt(diag(P));
    stdExtent = stdState(extendStart:end);
    phi = xUpd(3);
    [hull_x, hull_y] = pol2cart(model.gp.theta_test(:), xUpd(extendStart:end));
    [lower_hull_x, lower_hull_y] = pol2cart(model.gp.theta_test(:), xUpd(extendStart:end)-3*stdExtent);
    [upper_hull_x, upper_hull_y] = pol2cart(model.gp.theta_test(:), xUpd(extendStart:end)+3*stdExtent);
    hull_b = [hull_x' hull_x(1);
              hull_y' hull_y(1)];
    lower_hull_b = [lower_hull_x' lower_hull_x(1);
              lower_hull_y' lower_hull_y(1)];
    upper_hull_b = [upper_hull_x' upper_hull_x(1);
              upper_hull_y' upper_hull_y(1)];
    mean = rotate(hull_b, phi) + xUpd(1:2);
    lower = rotate(lower_hull_b, phi) + xUpd(1:2);
    upper = rotate(upper_hull_b, phi) + xUpd(1:2);

end

