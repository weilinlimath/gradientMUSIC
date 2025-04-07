%% Set-up
rng("default")
m         = 100;
supp      = 2*pi*[0.2,0.2+8/100,0.5];
sparsity  = length(supp);
% jitter    = (2*rand(1,sparsity)-1)/m; % add some small random jitter
% supp      = supp + jitter;
amp       = ones(sparsity,1); 
sigma     = 1;

%%
m1       = m-1;
Phi      = exp(1i*(-m1:m1)'*supp); % 2m-1 measurements
y        = Phi * amp; % clean measurement

% add noise
u = normrnd(0,sigma,2*m-1,1);
v = normrnd(0,sigma,2*m-1,1);
e = u+1i*v;

% Toeplitz matrix and SVD
T = Toeplitz_matrix(y);
[U,S,~] = svds(T,sparsity);
y_noisy = y + e; 
T_noisy = Toeplitz_matrix(y_noisy);
[U_noisy, ~, ~] = svds(T_noisy,sparsity);

% landscape associated with U
q_help1 = @(t) U'*exp(1i*(0:m-1)'*t)/sqrt(m);              
q = @(t) 1-(norm(q_help1(t),2))^2;       

% landscape associated with tilde U
q_help2   = @(t) U_noisy'*exp(1i*(0:m-1)'*t)/sqrt(m);      
q_noisy   = @(t) 1-(norm(q_help2(t),2))^2;  

theta = norm(U*U'-U_noisy*U_noisy',2);

%% plot
t = 2*pi*linspace(0,1,1e6);

q_samples = 0*t; 
q_noisy_samples = 0*t; 
for kk = 1:length(t)
    q_samples(kk) = q(t(kk));
    q_noisy_samples(kk) = q_noisy(t(kk));
end

figure;
plot(t,q_samples,LineWidth=2)
title('Graph of the landscape function associated with a $U$',FontSize=16,Interpreter='latex')
xlim([0,2*pi])
pbaspect([2*pi 1 1])
exportgraphics(gca,'landscapeplot1.pdf','BackgroundColor','none')

cmap = colororder();

figure;
plot(t,q_noisy_samples,LineWidth=2)
title('Graph of the landscape function associated with a $\widetilde U$',FontSize=16,Interpreter='latex')
xlim([0,2*pi])
pbaspect([2*pi 1 1])
exportgraphics(gca,'landscapeplot2.pdf','BackgroundColor','none')

figure;
plot(t,abs(q_samples-q_noisy_samples),LineWidth=2)
title('Graph of the absolute difference between landscape functions',FontSize=16,Interpreter='latex')
xlim([0,2*pi])
pbaspect([2*pi 1 1])
exportgraphics(gca,'landscapeplot3.pdf','BackgroundColor','none')

%% accepted and rejected elements

step = 1/(8*m); 
grid = 2*pi* ( 0:step:1 );

q_noisy_grid = 0*grid; 
for kk = 1:length(grid)
    q_noisy_grid(kk) = q_noisy(grid(kk));
end
idx = (q_noisy_grid <= 1/2); % accepted points on coarse grid

figure;
plot(t,q_noisy_samples,LineWidth=2,Color='k')
hold on 
plot(t,1/2*ones(1,length(t)),LineWidth=2,Color='k',LineStyle='--')
s1 = scatter(grid(idx),q_noisy_grid(idx),40,'filled');
s2 = scatter(grid(~idx),q_noisy_grid(~idx),40,'filled');
hold off
xlim([0,2*pi])
pbaspect([2*pi 1 1])
s1.MarkerEdgeColor = 'k';
s1.MarkerFaceColor = cmap(5,:);
s2.MarkerEdgeColor = 'k';
s2.MarkerFaceColor = cmap(2,:);
legend('$q_{\tilde U}$','$1/2$','Accepted','Rejected')
legend(Fontsize=32,Interpreter='latex',Location='southeast')
title('Thresholding the landscape function',FontSize=32,Interpreter='latex')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
exportgraphics(gca,'landscapeplot6.pdf','BackgroundColor','none')

%% zoomed in plot

a = supp(1)/(2*pi);

figure;
plot(t,q_samples,t,q_noisy_samples,LineWidth=2)
title('Zoomed-in graph of landscape functions',FontSize=16,Interpreter='latex')
exportgraphics(gca,'landscapeplot4.pdf','BackgroundColor','none')
xlim([min(2*pi*(a-1/m)),max(2*pi*(a+1/m))])
legend('$q_U$','$q_{\tilde U}$')
legend(Fontsize=16,Interpreter='latex',Location='southwest')
exportgraphics(gca,'landscapeplot4.pdf','BackgroundColor','none')

figure;
plot(t,q_noisy_samples,LineWidth=2,Color='k')
hold on 
plot(t,1/2*ones(1,length(t)),LineWidth=2,Color='k',LineStyle='--')
s1 = scatter(grid(idx),q_noisy_grid(idx),50,'filled');
s2 = scatter(grid(~idx),q_noisy_grid(~idx),50,'filled');
hold off
xlim([min(2*pi*(a-1/m)),max(2*pi*(a+1/m))])
s1.MarkerEdgeColor = 'k';
s1.MarkerFaceColor = cmap(5,:);
s2.MarkerEdgeColor = 'k';
s2.MarkerFaceColor = cmap(2,:);
legend('$q_{\tilde U}$','$1/2$','Accepted','Rejected')
legend(Fontsize=16,Interpreter='latex',Location='southeast')
title('Thresholding the landscape function',FontSize=16,Interpreter='latex')
exportgraphics(gca,'landscapeplot5.pdf','BackgroundColor','none')