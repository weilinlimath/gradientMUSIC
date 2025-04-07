%%
m = 100;                                % m_0
b = 4;                                  % beta
A = b/(b-1);                            % A(beta)
C0 = 1/6 - 2*A*E1(m,b,b) - 1/(6*m^2);   % C_0(m_0,beta)
C1 = 1/6 + 2*A*E1(m,b,b);               % C_1(m_0,beta)
theta = 0.01;                           % sin-theta dist

%% perturbation of x_j to tilde x_j

% calculate T_0, T_1, T_2, T_3 for this part
a = b-1/(2*pi);  % alpha for this part
T = Tconst(m,a,b);
T0 = T(1); T1 = T(2); T2 = T(3); T3 = T(4);

r = 7; 
if r*C0-1-(1/20 + T0*T3 + 3*T1*T2)*r^2*theta <= 0
    disp('invalid parameter choices')
end
% need this to be positive to use lemma

%% region of convexity

% calculate T_0, T_1 and T_2 for this step
a = b-1/6-r*theta/(2*pi);
T = Tconst(m,a,b);
T0 = T(1); T1 = T(2); T2 = T(3);

C2lower = (1-1/m^2)/6 - (pi/3)^2/30 - 2*T1^2 - 2*T0*T2 - theta - r*theta/10;
C2upper = 1/6 + 2*T1^2 + 2*T0*T2 + theta + r*theta/10;

%% lower and upper bound for first derivative

% calculate T_0 and T_1 for this step
a = b - 2/3 - r*theta/(2*pi);                       
T = Tconst(m,a,b);
T0 = T(1); T1 = T(2); 

% control derivative
lower = @(t) (1/3*(1-1/m^2) - t.^2/120) * sin(t/2);
C1lower = min( lower(pi/3), lower(4*pi/3) ) - 2*T0*T1 - theta - r*theta/6;
C1upper = 1/6 + 6/pi*T0*T1 + (theta + r*theta/6)*3/pi;

%% lower bound for q on far set

% calculate T_0 coonstant
T = Tconst(m,b/2,b);
T0 = T(1); 

% lower bound for far set
tau0 = 2/3-r*theta/(2*pi);
Cfar = 1 - max( 1/4, 1/(m^2*sin(pi*tau0/m)^2) ) - T0^2 - theta;

%% compute energy terms via numerical integration

function u = E0(m,a,b)
g = @(x) 1./(m*sin(x/2));
h = @(x) g(x).^2;
u = integral(h,2*pi*(a+b/2)/m,pi);
u = 2*u*m/(2*pi*b)+2*h(2*pi*a/m);
end

function u = E1(m,a,b)
g = @(x) 1./(m*sin(x/2));
h = @(x) (g(x)/2+g(x).^2/2).^2;
u = integral(h,2*pi*(a+b/2)/m,pi);
u = 2*u*m/(2*pi*b)+2*h(2*pi*a/m);
end

function u = E2(m,a,b)
g = @(x) 1./(m*sin(x/2));
h = @(x) (g(x)/4+g(x).^2/2+g(x).^3/2).^2;
u = integral(h,2*pi*(a+b/2)/m,pi);
u = 2*u*m/(2*pi*b)+2*h(2*pi*a/m);
end

function u = E3(m,a,b)
g = @(x) 1./(m*sin(x/2));
h = @(x) (g(x)/8+3*g(x).^2/8+6*g(x).^3/8+10*g(x).^4/8).^2;
u = integral(h,2*pi*(a+b/2)/m,pi);
u = 2*u*m/(2*pi*b)+2*h(2*pi*a/m);
end

function u = Tconst(m,a,b)
u = zeros(1,4);
A = b/(b-1);
u(1) = sqrt(A*E0(m,a,b)) + 2*sin(asin(sqrt(A*E0(m,b,b)))/2);
u(2) = sqrt(A*E1(m,a,b)) + 2*sin(asin(sqrt(A*E0(m,b,b)/12))/2);
u(3) = sqrt(A*E2(m,a,b)) + 2*sin(asin(sqrt(A*E0(m,b,b)/80))/2);
u(4) = sqrt(A*E3(m,a,b)) + 2*sin(asin(sqrt(A*E0(m,b,b)/448))/2);
end