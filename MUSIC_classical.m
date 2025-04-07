function supp_classical = MUSIC_classical(U_noisy, width)

[L, sparsity] = size(U_noisy); 
Grid        = (0 : width : 2*pi-width)';
GridLen     = length(Grid);

%funphi       = @(t) exp(1i*(0:L-1)'*t);
funs_noisy   = @(t) U_noisy'*exp(1i*(0:L-1)'*t)/sqrt(L);                          % function r(t) 
funsd_noisy  = @(t) U_noisy'*((1i*(0:L-1)').*exp(1i*(0:L-1)'*t))/sqrt(L);         % function r'(t)
funsd2_noisy = @(t) U_noisy'*((-((0:L-1).^2)').*exp(1i*(0:L-1)'*t))/sqrt(L);         % function r''(t)
funq_noisy   = @(t) 1-(norm(funs_noisy(t),2))^2;  %norm(NoiseSpace_noisy'*exp(1i*(0:L-1)'*t),2)^2/L;                      % function q(t)
funqd_noisy  = @(t) -2*real((funs_noisy(t))'*funsd_noisy(t)); %sum(2*real(funr_noisy(t).*conj(funrd_noisy(t))));         % function q'(t)


%% noise subspace correlation function
ImagingFun_noisy = zeros(GridLen,1);
%RFun_noisy  = 1 - vecnorm(SignalSpace_noisy'*exp(1i*(0:L-1)'*Grid)/sqrt(L)).^2;
parfor k = 1 : GridLen
    ImagingFun_noisy(k) = 1/funq_noisy(Grid(k));       %(norm(NoiseSpace_noisy'*funphi_noisy(Grid(k)),2))^2/L;
end


%% Step 2 - Compute local maxima
LocalCompareMatrix = zeros(GridLen,2);
LocalCompareMatrix(2:end,1)   = ImagingFun_noisy(2:end)-ImagingFun_noisy(1:end-1);
LocalCompareMatrix(1,1)       = ImagingFun_noisy(1)-ImagingFun_noisy(end);
LocalCompareMatrix(1:end-1,2) = ImagingFun_noisy(1:end-1)-ImagingFun_noisy(2:end);
LocalCompareMatrix(end,2)     = ImagingFun_noisy(end)-ImagingFun_noisy(1);
LocalCompareValue             = min(LocalCompareMatrix,[],2);
LocalMaxIndex                 = find(LocalCompareValue>0);
LocalMaxImagingFun_noisy      = ImagingFun_noisy(LocalMaxIndex);
[~, LocalMaxSLargestInd]      = maxk(LocalMaxImagingFun_noisy,sparsity);
supp_classical = Grid(LocalMaxIndex(LocalMaxSLargestInd));





