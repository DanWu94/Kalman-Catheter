function R = R_euler(a,b,c)
%Calculate Euler Angle
R = [cos(a)*cos(c)-sin(a)*sin(c)*cos(b),sin(a)*cos(b)-sin(a)*sin(c)*cos(b),sin(c)*sin(b);
    -cos(a)*sin(c)-sin(a)*cos(c)*cos(b),-sin(a)*sin(c)+cos(a)*cos(b),cos(c)*sin(b);
    sin(a)*sin(b),-cos(a)*sin(b),cos(b)];
end

