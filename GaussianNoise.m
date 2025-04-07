%% Error versus number of measurements

%% Start parallel pool if not initialized
if isempty(gcp('nocreate'))
    parpool("Processes")
end

%% Set-up
sparsity  = 5;         % sparsity or number of frequencies
amp       = ones(sparsity,1); %.*exp(2*pi*1i*rand(sparsity,1));
amin      = min(abs(amp));
Marray    = round(10.^(2:.5:4));
MLen      = length(Marray);
Ntrial    = 50;
sep       = 8*pi/min(Marray);
jitter    = 0.5*(2*rand(1,sparsity)-1); % add some small random jitter
supp      = 0.2 + 2*sep*((1:sparsity)+jitter)';
%supp      = sort(supp,'ascend');

% noise parameters
sigma    = 0.5;
Rgrowth  = [-1/4,0,1/4];
Ngrowth   = length(Rgrowth);  

%% Initialize matrices

XErrorGDMUSIC   = zeros(MLen,Ngrowth,Ntrial);
AErrorGDMUSIC  = zeros(MLen,Ngrowth,Ntrial);

%% main loop for experiment 
for iterm = 1:MLen
    m = Marray(iterm)  % bandwidth
    m1 = m-1; 
    Grid     = 0 : (0.5/m) : 2*pi;     % gradient MUSIC grid
    Phi      = exp(1i*(-m1:m1)'*supp'); % 2m-1 measurements
    y        = Phi * amp; % clean measurement
    
    for iterr = 1:Ngrowth
        
        % grow/decay rate of noise
        r = Rgrowth(iterr); 

        for itern = 1 : Ntrial
            
            % add noise
            u = ((1+abs(-m1:m1)).^r)'.*normrnd(0,sigma,2*m-1,1);
            v = ((1+abs(-m1:m1)).^r)'.*normrnd(0,sigma,2*m-1,1);
            e = u+1i*v; 
            y_noisy  = y + e;
    
            % Toeplitz matrix and SVD
            T_noisy = Toeplitz_matrix(y_noisy);
            [U_noisy, ~, ~] = svds(T_noisy,sparsity);
    
            % Gradient MUSIC parameters
            threalpha   = 0.529;       % threshold
            GDstepsize  = 5/m^2; 
            GDtol       = 10^(-2);
            GDiter      = [31,300];
    
            % Gradient MUSIC
            supp_gdmusic = MUSIC_gradient(U_noisy,Grid,threalpha,GDstepsize,GDtol,GDiter);
            supp_gdmusic = sort(supp_gdmusic, 'ascend');
            amp_gdmusic = amp_quadratic(T_noisy, supp_gdmusic);
            error_gdmusic = dist_parameter(supp, amp, supp_gdmusic, amp_gdmusic);
    
            % update matrices
            XErrorGDMUSIC(iterm,iterr,itern) = error_gdmusic(1);
            AErrorGDMUSIC(iterm,iterr,itern)= error_gdmusic(2);
    
        end
    end
end


%% compute the 90% quantile over trials

XErrorMGDMUSIC = quantile(XErrorGDMUSIC,0.9,3);
AErrorMGDMUSIC= quantile(AErrorGDMUSIC,0.9,3);

%% support error plot
figure;
loglog(Marray,XErrorMGDMUSIC,Marker='*',LineWidth=2,MarkerSize=8)

legendhelp = cell(1,Ngrowth);
slope_gd_music = zeros(1,Ngrowth);
for iterr = 1:Ngrowth
    poly_music = polyfit(log10(Marray),log10(XErrorMGDMUSIC(:,iterr)),1);
    slope_gd_music(iterr) = poly_music(1);
    legendhelp{iterr} = ['Gradient-MUSIC slope $\approx$ ', num2str(slope_gd_music(iterr),'%6.3f'), ', $r=' num2str(Rgrowth(iterr)) '$']; 
end

legend(legendhelp,FontSize=16,Interpreter='latex',Location='southwest')
title('Support error versus $m$',FontSize=16,Interpreter='latex')
xlabel('Bandwidth $m$',FontSize=16,Interpreter='latex')
ylabel('Support Error',FontSize=16,Interpreter='latex')
set(gca,'FontSize', 16)

exportgraphics(gca,'GaussianNoiseSupport.pdf','BackgroundColor','none')

%% amplitude error plot
figure;
loglog(Marray,AErrorMGDMUSIC,Marker='*',LineWidth=2,MarkerSize=8)

legendhelp = cell(1,Ngrowth);
slope_gd_music = zeros(1,Ngrowth);
for iterr = 1:Ngrowth
    poly_music = polyfit(log10(Marray),log10(AErrorMGDMUSIC(:,iterr)),1);
    slope_gd_music(iterr) = poly_music(1);
    legendhelp{iterr} = ['Gradient-MUSIC slope $\approx$ ', num2str(slope_gd_music(iterr),'%6.3f'), ', $r=' num2str(Rgrowth(iterr)) '$']; 
end

legend(legendhelp,FontSize=16,Interpreter='latex',Location='southwest')
title('Amplitude error versus $m$',FontSize=16,Interpreter='latex')
xlabel('$m$',FontSize=16,Interpreter='latex')
ylabel('Amplitude error',FontSize=16,Interpreter='latex')
set(gca,'FontSize',16)

exportgraphics(gca,'GaussianNoiseAmplitude.pdf','BackgroundColor','none')

%%
filename = ['GaussianNoiseTrial' num2str(Ntrial) '.mat' ];
save(filename)
