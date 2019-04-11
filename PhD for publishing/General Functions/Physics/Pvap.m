function [ P ] = Pvap( T )
%This function computes the saturated vapor pressure corresponding to a
%given temperature
% P (Pa)
% T (K)

P = 3385.5*exp(-8.0929+0.97608*((T+42.607-273.15)^(0.5)));

end

