function a_rec = amp_quadratic(T,supp)
% given (noisy) Toeplitz matrix and estimated support, compute estimated
% amplitudes using quadratic method

[m,n] = size(T); 
if m~=n 
    error('Toeplitz matrix not square')
end
supp = supp(:); 

m = m-1;
Phi = exp(1i*(-m/2:m/2)'*supp');
Phi_inv = pinv(Phi); 
a_rec = diag(Phi_inv*T*Phi_inv');  