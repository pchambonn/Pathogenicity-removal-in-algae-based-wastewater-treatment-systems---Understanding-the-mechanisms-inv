function [ T ] = Tvap( P )
%This function computes the temperature of the air from the value of
%saturated water vapor pressure

%% Inputs:

% P: saturated water vapor pressure (Pa)

%% Outputs

% T: air temperature (°C)

%% Calculation

T = ((1/0.97608)*(8.0929 + log(P/3385.5)))^(2) - 42.607;

end

