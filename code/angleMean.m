function mean = angleMean(angles,w)
    wsum = sum(w);
    c = sum(w'.*cos(angles))./wsum;
    s = sum(w'.*sin(angles))./wsum;
    mean = ssa(atan2(s,c));
end