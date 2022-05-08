function [theta,phi,r] = cart_to_sph(x,y,z)
 hypotxy = hypot(x,y);
 r = hypot(hypotxy,z);
 theta = atan2(hypotxy,z);
 phi = atan2(y,x);
 end

