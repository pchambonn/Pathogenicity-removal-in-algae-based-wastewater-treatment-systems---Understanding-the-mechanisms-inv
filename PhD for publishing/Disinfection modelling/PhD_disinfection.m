function [ Kd ] = PhD_dinsinfection( t , Tp , pH , G0 , sigma_a , d , alpha_dis )

%% This functions computes the disinfection rate (first order)according to my PhD findings.

%% INPUTS:

% t: period of time during which the decay rate is calculated (s) 
% Tp: pond temperature (K)
% pH: pond pH
% G0: Solar irradiance (W/m2)
% sigma_a: algae abosrbance (g.m-2)
% d: pond depth (m)
% alpha_dis: empiric linear law of disinfection against direct sunlight

%% OUTPUTS:

% Kd: decay rate (d-1)

k_dark = 10.4*1.14^(Tp - 273.15 - 20);

k_pH = 3.1*10^4*1.14^(Tp - 273.15 - 20)*10^(pH - 14);

k_s = alpha_dis*G0/sigma_a/d*(1 - exp(-sigma_a*d));

Kd = k_dark + k_pH + k_s;


end

