function IOU = IntersectOverUnion(hull1,hull2)
    %calculates the IOU measure of two closed 2D shapes (Probably a nicer
    %way to handle this but quick and easy for now)
    poly1 = polyshape(hull1);
    poly2 = polyshape(hull2);
    IOU = area(intersect(poly1,poly2))/area(union(poly1,poly2));
end