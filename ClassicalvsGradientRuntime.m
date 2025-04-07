%% Error versus number of measurements

%% Start parallel pool if not initialized
if isempty(gcp('nocreate'))
    parpool("Processes")
end

%% Set-up
sparsity  = 5;         % sparsity or number of frequencies
amp       = ones(sparsity,1); %.*exp(2*pi*1i*rand(sparsity,1));
amin      = min(abs(amp));
Marray    = 10^3;
MLen      = length(Marray);
Ntrial    = 10;
sep       = 8*pi/min(Marray);
jitter    = 0.5*(2*rand(1,sparsity)-1); % add some small random jitter
supp      = 0.2 + 2*sep*((1:sparsity)+jitter)';

% noise parameters
sigma    = 0.1;
r        = 0; % white noise

%% Initialize matrices

TErrorMNGDMUSIC = zeros(MLen,Ntrial);
TErrorMNMUSIC   = zeros(MLen,Ntrial);
TSVDMN          = zeros(MLen,Ntrial);

%% main loop for experiment 
for iterm = 1:MLen
    m = Marray(iterm);  % bandwidth
    m1 = m-1; 
    h = 0.1*sigma/m^(3/2-r); % classical MUSIC grid width
    Grid     = 0 : (0.5/m) : 2*pi;     % gradient MUSIC grid
    Phi      = exp(1i*(-m1:m1)'*supp'); % 2m-1 measurements
    y        = Phi * amp; % clean measurement

    % add noise
    for itern = 1 : Ntrial
        u = ((1+abs(-m1:m1)).^r)'.*normrnd(0,sigma,2*m-1,1);
        v = ((1+abs(-m1:m1)).^r)'.*normrnd(0,sigma,2*m-1,1);
        e = u+1i*v; % blue or red noise
        y_noisy  = y + e;

        % Toeplitz matrix and SVD
        tic;
        T_noisy = Toeplitz_matrix(y_noisy);
        [U_noisy, ~, ~] = svds(T_noisy,sparsity);
        time_svd = toc;
        TSVDMN(iterm,itern) = time_svd;

        % Gradient MUSIC parameters
        threalpha   = 0.529;       % threshold
        GDstepsize  = 5/m^2; 
        GDtol       = 10^(-2);     % GDtol corresponds to ep m 
        GDiter      = [31,300];

        % Gradient MUSIC
        tic; 
        supp_gdmusic = MUSIC_gradient(U_noisy,Grid,threalpha,GDstepsize,GDtol,GDiter);
        time_gdmusic = toc;
        supp_gdmusic = sort(supp_gdmusic, 'ascend');
        amp_gdmusic = amp_quadratic(T_noisy, supp_gdmusic);
        error_gdmusic = dist_parameter(supp, amp, supp_gdmusic, amp_gdmusic);

        % Classical MUSIC
        tic
        supp_music  = MUSIC_classical(U_noisy, h);
        time_music  = toc;
        supp_music  = sort(supp_music,'ascend');
        amp_music = amp_quadratic(T_noisy, supp_music);
        error_music = dist_parameter(supp, amp, supp_music, amp_music);

        % update matrices
        TErrorMNGDMUSIC(iterm,itern) = time_gdmusic;
        TErrorMNMUSIC(iterm,itern) = time_music;
        TSVDMN(iterm,itern) = time_svd;
    end
end


%% compute the 90% quantile over trials

TErrorMGDMUSIC = max(TErrorMNGDMUSIC,[],2);
TErrorMMUSIC = max(TErrorMNMUSIC,[],2);
TSVDM = max(TSVDMN,[],2);

%%
% filename = ['RuntimeTrial' num2str(Ntrial) '.mat' ];
% save(filename)
