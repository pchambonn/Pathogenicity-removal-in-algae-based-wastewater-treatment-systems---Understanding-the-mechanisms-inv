function [ hconv ] = conv_coeff( Nu_a,alpha_a,L,Nu,lambda_a )
% This function computes the convection coefficient at the surface of an
% opaque pond of a given geometry

% Inputs
% lambda_a : air thermal conductivity (W/m/K)
% alpha_a : air thermal diffusivity (m²/s)
% Nu_a : air kinematic viscosity (m²/s)
% Nu : wind velocity (m/s)
% L : characteristic pond length (m)

% Prandtl
Pr = Nu_a/alpha_a ;

% Reynolds
Re = L*Nu/Nu_a ;

if Re > 5*10^5
    A = 0.035 ;
    beta = 0.8 ;   
    delta = 1/3 ;
    hconv = lambda_a/L*A*(Re^beta)*(Pr^delta) ;
else
    if Re > 3*10^5
        h_conv_1 = lambda_a/L*0.035*(Re^0.8)*(Pr^(1/3));
        h_conv_2 = lambda_a/L*0.628*(Re^0.5)*(Pr^(1/3));
        hconv = h_conv_2 + (h_conv_1 - h_conv_2)*(Re - 3*10^5)/(5*10^5 - 3*10^5);
    else
        A = 0.628 ;
        beta = 0.5 ;   
        delta = 1/3 ;
        hconv = lambda_a/L*A*(Re^beta)*(Pr^delta) ;
    end
end

end

