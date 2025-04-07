function T = Toeplitz_matrix(y,m)
% creates a square Toeplitz matrix of y with m+1 rows

if nargin < 2
    m = floor((length(y)-1)/2);
end

if 2*m+1 > length(y)
    disp('invalid parameter choice')
end

y = y(1:2*m+1); 
T = zeros(m+1,m+1); 

for kk = 1:m+1
    T(:,kk) = y(m+1-kk+1:2*m+1-kk+1);
end
