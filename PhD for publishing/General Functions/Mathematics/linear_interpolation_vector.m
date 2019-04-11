function [ B ] = linear_interpolation_vector( A , p )
%This functions aims to create from a vector of size n a new vector of size
%p*n of which coordinates are the linear interpolation of the starting
%vector's coordinates.

n = length(A);
B = zeros(1,p*(n-1)+1);

for k =1:n-1
    for i = 1:p
        B(p*(k-1) + i) = A(k) + (i-1)*(A(k+1) - A(k))/p ;
    end
end

B(1,p*(n-1)+1) = A(n);

end

