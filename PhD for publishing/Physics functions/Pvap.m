function [ P ] = Pvap( T )
%This function computes the saturated water vapor pressure corresponding at a
%given temperature at Earth surface

%% Inputs

% T: temperature (K)

%% Outputs

% P: water apor pressure (Pa)

%% Calculation

P = 3385.5*exp(-8.0929+0.97608*((T+42.607-273.15)^(0.5)));

end

