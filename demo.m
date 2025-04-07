%% demo for code

%% Start parallel pool if not initialized
if isempty(gcp('nocreate'))
    parpool("Processes")
end

%% Create parameters and signal

sparsity  = 5; 
m = 100; % bandwidth
m1 = m-1; 
amp = ones(sparsity,1); % amplitudes
sep = 8*pi/m;
supp = (1:sparsity)'*sep;
supp = sort(supp,'ascend'); % frequencies satisfying separation condition
Phi = exp(1i*(-m1:m1)'*supp'); % 2m+1 measurements
y = Phi * amp; % clean measurement

% clean Toeplitz matrix
T = Toeplitz_matrix(y); 
[U,S,~] = svds(T, sparsity); % empirical signal space

% Create some noise
sigma = 0.1;
r = 0;
e = ((1+abs(-m1:m1)).^r)'.*normrnd(0,sigma,2*m-1,1);
y_noisy  = y + e; % noisy data

%% Gradient MUSIC full pipeline

[supp_gdmusic2, amp_gdmusic2] = MUSIC_gradient_full(y_noisy);

%% Demos for approximation versions of algorithms

% Toeplitz matrix and SVD
T_noisy = Toeplitz_matrix(y_noisy);
[U_noisy,~,~] = svds(T_noisy, sparsity); % empirical signal space

% verify factorization of  T
Phi2 = exp(1i*(0:m1)'*supp');
norm(T-Phi2*diag(amp)*Phi2');

% compute vartheta and rho 
theta = norm(U*U'-U_noisy*U_noisy',2);
rho = 2*norm(Toeplitz_matrix(e),2)/(S(sparsity,sparsity));


%% classical MUSIC 

% parameters
mesh = 0.01/m^(3/2-r); % width of a uniform grid
supp_music = MUSIC_classical(U_noisy,mesh);
supp_music = sort(supp_music,'ascend');
amp_music = amp_quadratic(T_noisy,supp_music);
error_music = dist_parameter(supp,amp,supp_music,amp_music);

%% Gradient MUSIC core approximation

% parameters
mesh        = 1/(2*m); % width of uniform grid
Grid        = 0 : mesh : 2*pi;
threalpha   = 0.529;       % threshold
GDtol       = 10^(-2);
GDiter      = [31, 300];
GDstepsize  = 5/m^2; 

supp_gdmusic = MUSIC_gradient(U_noisy,Grid,threalpha,GDstepsize,GDtol,GDiter); 
supp_gdmusic  = sort(supp_gdmusic,'ascend');
amp_gdmusic = amp_quadratic(T_noisy,supp_gdmusic);
error_gdmusic = dist_parameter(supp,amp,supp_gdmusic,amp_gdmusic);