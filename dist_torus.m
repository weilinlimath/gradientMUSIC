function d = dist_torus(x,y)
% computes the distance between arrays x and y
% the j-th entry of d is the distance between x_j and y_j

x = mod(x,2*pi);
y = mod(y,2*pi);
d = min(min(abs(x-y),abs(x-y+2*pi)),abs(x-y-2*pi));