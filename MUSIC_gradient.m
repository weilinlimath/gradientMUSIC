function supp_gdmusic = MUSIC_gradient(U_noisy,Grid,threalpha,GDstepsize,GDtol,GDiter)

GDitermin = GDiter(1);
GDitermax = GDiter(2)-GDiter(1);
m = size(U_noisy,1);
GridLen = length(Grid);

funs_noisy   = @(t) U_noisy'*exp(1i*(0:m-1)'*t)/sqrt(m);                          % function r(t)
funsd_noisy  = @(t) U_noisy'*((1i*(0:m-1)').*exp(1i*(0:m-1)'*t))/sqrt(m);         % function r'(t)
funq_noisy   = @(t) 1-(norm(funs_noisy(t),2))^2;                                  % function q(t)
funqd_noisy  = @(t) -2*real((funs_noisy(t))'*funsd_noisy(t));                     % function q'(t)

% noise subspace
Grid_q_noisy = zeros(GridLen,1);
for k = 1 : GridLen
    Grid_q_noisy(k) = funq_noisy(Grid(k));
end

%% Step 1 - Threshold the landscape function and find representatives

ThreIndex = double(Grid_q_noisy<threalpha);

% find endpoints of the clusters
ThreIndexR = ThreIndex-[ThreIndex(2:end);ThreIndex(1)];
Rcluster = Grid(ThreIndexR==1);
ThreIndexL = ThreIndex-[ThreIndex(end);ThreIndex(1:end-1)];
Lcluster = Grid(ThreIndexL==1);

% deal with potential wrap around effect
if dist_torus(Rcluster(1),Lcluster(1))>2*pi/m && dist_torus(Rcluster(1),Lcluster(end))<2*pi/m
    Lcluster = [Lcluster(end), Lcluster(1:end-1)];
end

NumInt = length(Rcluster);
GDInitial = mod(Lcluster + dist_torus(Lcluster,Rcluster)/2, 2*pi);

% ThreGridGap      = ThreGrid(2:end) - ThreGrid(1:end-1);
% ThreGridGapInd   = find(ThreGridGap>=1.5*mesh);
% NumInt           = 0;
% if isempty(ThreGrid) == 0
%     ThreGridInt      = [];
%     ThreGridInt(1,1) = ThreGrid(1);
%     for k = 1:length(ThreGridGapInd)
%         ThreGridInt(k,2)   = ThreGrid(ThreGridGapInd(k));
%         ThreGridInt(k+1,1) = ThreGrid(ThreGridGapInd(k)+1);
%     end
%     ThreGridInt(k+1,2) = ThreGrid(end);
%     NumInt             = length(ThreGridGapInd)+1;
%     clear ThreGridGap
%     % combine the first interval and the last interval if they are connected
%     ifwrap = 0;
%     if abs(ThreGridInt(1,1)-ThreGridInt(NumInt,2)+2*pi)<=2*mesh
%         ThreGridInt(1,1)      = ThreGridInt(NumInt,2);
%         ThreGridInt(NumInt,:) = [];
%         NumInt = NumInt - 1;
%         ifwrap = 1;
%     end
% end

%% Step 2 - Gradient Descent
supp_gdmusic = zeros(NumInt,1);

parfor k = 1 : NumInt
    iter = 0;
    % Gradient descent for the k-th interval
    xiter   = GDInitial(k);
    xgrad   = funqd_noisy(xiter);
    for jj = 1:GDitermin
        xnext = mod(xiter - GDstepsize * xgrad, 2*pi);
        xiter = xnext;
        xgrad   = funqd_noisy(xiter);
    end
    while (abs(xgrad) > GDtol) && (iter <= GDitermax)
        iter  = iter + 1;
        xnext = mod(xiter - GDstepsize * xgrad, 2*pi);
        xiter = xnext;
        xgrad   = funqd_noisy(xiter);
    end
    supp_gdmusic(k)  = xiter;
end
