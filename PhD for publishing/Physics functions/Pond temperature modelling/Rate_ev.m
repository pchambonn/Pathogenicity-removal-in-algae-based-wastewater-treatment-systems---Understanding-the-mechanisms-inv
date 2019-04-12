function [ me ] = Rate_ev( L,Dw_a,Nu_a,Nu,Tp,RH,R,Ta,Mw )

% This function calculates the rate of evaporation of any given opaque pond


%% Inputs
% L : characteristic pond length (m)
% Dw_a : mass difusion coefficiet of water vapor in air (m²/s)
% Nu_a : air kinematic viscosity (m²/s)
% Nu : wind velocity (m/s)
% Tp : pond temperature (K)
% RH : relative humidity of the air above the pond surface
% R : ideal gas constant (Pa.m3/mol/K)
% Ta : air temperature (K)
% Mw : molecular weight of water (kg/mol)

%% Outputs

% me: rate of evaporation (kg/m2/s)

%% Calculation

% Schmidt number:
Sch = Nu_a/Dw_a ;

% Reynolds number :
Re = L*Nu/Nu_a ;

if Re > 5*10^5
    A = 0.035 ;
    beta = 0.8 ;   
    delta = 1/3 ;
    K = Dw_a/L*A*(Re^beta)*(Sch^delta) ;
else 
    if Re > 3*10^5
        K_1 = Dw_a/L*0.035*(Re^0.8)*(Sch^(1/3));
        K_2 = Dw_a/L*0.628*(Re^0.5)*(Sch^(1/3));
        K = K_2 + (K_1 - K_2)/(5*10^5 - 3*10^5)*(Re - 3*10^5);
    else
        A = 0.628 ;
        beta = 0.5 ;   
        delta = 1/3 ;
        K = Dw_a/L*A*(Re^beta)*(Sch^delta) ;
    end
end

me = K*(Pvap(Tp)/Tp - RH*Pvap(Ta)/Ta)*Mw/R;

end

