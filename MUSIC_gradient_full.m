function [supp_gdmusic, amp_gdmusic] = MUSIC_gradient_full(y_noisy)
% Full pipeline version
% all parameters are chosen accorinding to the theorems

%% Step 1: sparsity detection and Toeplitz estimator

% parameters
gamma = 0.0525;

T_noisy = Toeplitz_matrix(y_noisy);
[U,S,~] = svd(T_noisy);
S = diag(S);
Sratio = S/S(1);
sparsity = sum(double(Sratio>=gamma));
U_noisy = U(:,1:sparsity);
m = size(U_noisy,1);

%% Step 2: use core approximation version of gradient-MUSIC

% Gradient MUSIC parameters
alpha       = 0.529;       % threshold
GDstepsize  = 5/m^2;
GDtol       = 10^(-2);     % GDtol corresponds to ep m
GDiter      = [31,300];
Grid        = 0 : (0.5/m) : 2*pi;     % gradient MUSIC grid

supp_gdmusic = MUSIC_gradient(U_noisy,Grid,alpha,GDstepsize,GDtol,GDiter);
supp_gdmusic = sort(supp_gdmusic, 'ascend');

%% Step 3: use quadratic method to estimate amplitudes

amp_gdmusic = amp_quadratic(T_noisy, supp_gdmusic);

