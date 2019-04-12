function [ CO2_sat ] = CO2_solubility( T )

% This functions aims at calculating the solubility of carbon dioxide in
% water under atmospheric pressure according to temperature T (K)
% The data used is from Metcalf and Eddy for the Henry constant law of
% temperature (using a linear interpolation), assuming 400 ppm of CO2 in
% the atmosphere and a constant air pressure of 1 atm.
% H (Henry's constant) is expressed in atm.
% p_CO2 (mole fraction)

%% INPUT:

% T : water temperature (K)

%% OUTPUT

% CO2_sat: carbon dioxide saturation concentration (g C-CO2/m3)


%% Calculations

H = 10^(6.606 - 1012.40/T);
CO2_sat = 1/(1.50*10^(-6))*1/H*400*10^(-6);

end

