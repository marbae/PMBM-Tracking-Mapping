function newPositions = rotate(positions, angle)
    % rotates positions by an angle
    % positions, newPositions : 2 x Nt
    % angle: 1
    c = cos(angle);
    s = sin(angle);
    R = [c -s; s c];
    newPositions = R*positions;
end
