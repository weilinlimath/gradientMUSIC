function d = dist_parameter(supp1,amp1,supp2,amp2)
% computes the error between supports after they have been best matched
% and amplitude error according to this best matching. Error is taken in
% the sup norm.

% all inputs must be of the same length
if length(supp1) == length(supp2) && length(amp1) == length(amp2) ...
        && length(supp1)==length(amp1)
    sparsity = length(supp1);

    % compute all permutations
    Permutations = perms(1:sparsity);
    Nperm        = size(Permutations,1);
    DistPerm     = zeros(Nperm,1);
    for k = 1:Nperm
        supp2permuted = supp2(Permutations(k,:));
        DistPerm(k) = max(dist_torus(supp1,supp2permuted));
    end

    [supp_error, index] = min(DistPerm);
    bestP = Permutations(index,:);
    amp_error = norm(amp1-amp2(bestP),inf);
    d = [supp_error, amp_error];

else
    disp('recovered sparsity does not match true sparsity')
    d = [pi, inf]; 
end
