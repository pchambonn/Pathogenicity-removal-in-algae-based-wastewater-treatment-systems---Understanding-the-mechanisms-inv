function [ B ] = refinement_linear_interpol( A, n_ref, i )
% This function aims at refining calculations by expanding a given table A
% with n_ref values between the i-eme and i+1-eme coordinates. The new
% values are the n_ref linear interpolation of A(i);A(i+1).
% NB: will only work with a vector i.e. size ~ (1,n)

B = A;

for k = 1:n_ref - 1
    B = [B(1:i + k - 1),A(i) + k/n_ref*(A(i+1) - A(i))];
end

B = [B, A(i + 1:end)];

end

