function [ T ] = Tvap( P )
%This function computes the temperature of the air from the saturated
%vapor pressur
% P (Pa)
% T (*C)

T = ((1/0.97608)*(8.0929 + log(P/3385.5)))^(2) - 42.607;

end

