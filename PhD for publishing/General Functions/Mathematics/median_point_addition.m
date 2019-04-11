function [ B ] = median_point_addition( A , k )
% This functions aims at inserting a data point in a vector at a spot n that is the median point between the point k and k + 1. The size of the vector in incremented by one in the process. 

%% INPUTS

% A: vector (size(1,n) or size(n,1) regardless
% k data point at which a point should be added

if size(A,1) == 1
    B = [A(1:k) , (A(k) + A(k+1))/2 , A(k+1:end)];
end
if size(A,2) == 1
    B = [A(1:k) ; (A(k) + A(k+1))/2 ; A(k+1:end)];
else


end

