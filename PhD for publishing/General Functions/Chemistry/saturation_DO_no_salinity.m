function [ DO ] = saturation_DO_no_salinity( T )

%% This function aims at calculating the dissolved oxygen in a water body with no salinity

% Inputs

% T: temperature of water (K)

% Outputs

% DO: dissolved oxygen at saturation in the water body (mg/L)

%% Function execution

DO = 1.42905*exp(-173.4292 + 249.6339*(100/T) + 143.3483*log(T/100) - 21.8492*(T/100));

end

