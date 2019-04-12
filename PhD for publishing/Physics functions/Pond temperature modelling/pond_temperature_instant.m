function [ Tp,  me, B, Qra_p, Qra_s, Qev, Qconv,  Qcond, Qi, Qr ] = pond_temperature_instant(Tp,Ta,Hs_dir,Hs_dif,RH,qr,Dw_a,Nu,T_soil,Ti,qi,s_t,V,S,L,rho_w,Cp_w,eps_w,sigma_stefan,fa,eps_a,Nu_a,R,Mw,Lw,lambda_a,alpha_a,alpha_s,Ts_ref,ks)
% This function aims at calculating the temperature of a pond after a given
% step of time s_t starting from known conditions


%% Inputs

% Tp: intial temperature (K)
% Ta: air temperature (K)
% Hs_dir: direct sunlight intensity (W.m-2)
% Hs_dif: diffuse sunlight intensity (W.m-2)
% RH: relative humidity (-)
% qr: rainfall (mm)
% Dw_a : mass difusion coefficiet of water vapor in air (m²/s)
% Nu: wind speed (m.s-1)
% T_soil: matrix of soil temeprature profile (K)
% Ti : temperature of inlet (°C)
% qi : inflow rate (m3/s)
% s_t: step of time for the calculation of next pond temperature (s)
% V: pond volume (m3)
% S : pond surface (m²)
% L : characteristic pond length (m)
% rho_w : density of pond water (kg/m3)
% Cp_w : speciffic heat capacity of pond water (J/kg/K)
% eps_w : water emissivity
% sigma_stefan : Stefan-Boltzmann constant (W/m²/K) 
% fa : algal absorption fraction (%)
% eps_a : air emissivity
% Nu_a : air kinematic viscosity (m²/s)
% R : ideal gas constant (Pa.m3/mol/K)
% Mw : molecular weight of water (kg/mol)
% Lw : water latent heat (J/kg)
% lambda_a : air thermal conductivity (W/m/k)
% alpha_a : air thermal diffusivity (m²/s)
% alpha_s: soil thermal diffusivity (m2/s)
% Ts_ref: soil reference temperature (K)
% ks : soil thermal conductivity (W/m/K)

%% Outputs

% Tp: pond calculated temperature (K)
% me: rate of evaporation (kg/m2/s)
% Qra_p : heat loss through pond radiation (W)
% Qra_s : heat gain through solar radiation (W)
% Qra_a : heat gain through air radiation (W)
% Qev : heat loss through evaporation (W)
% me : rate of evaporation (kg/s/m²)
% Qconv : heat exchange through convection (W)
% Qcond: heat exchange through conduction (W)
% Qi:  heat exchange through inlet/outlet (W)
 % Qr: heat exchange through rainfall (W)


Qra_p = - eps_w*sigma_stefan*Tp^4*S ;

% Solar radiation

Qra_s = (1-fa)*(Hs_dir+Hs_dif)*S ;

% Air radiation

Qra_a = eps_w*eps_a*sigma_stefan*Ta^4*S ;

% Evaporation

me = Rate_ev(L,Dw_a,Nu_a,Nu,Tp,RH,R,Ta,Mw);

Qev = - me*Lw*S ;

% Convection

hconv = conv_coeff(Nu_a,alpha_a,L,Nu,lambda_a); 
Qconv = hconv*(Ta-Tp)*S ;

% Conduction

ls_ref = 4400*sqrt(alpha_s) ;

%Résolution de B pour la colonne i+1
p = size(T_soil,1);
B = zeros(p,1);
A = alpha_s*s_t/((ls_ref/p)^2);

B(1) = Tp;
B(p) = Ts_ref;

for k = 2:p-1 
    B(k) = A*(T_soil(k+1)-2*T_soil(k) + T_soil(k-1)) + T_soil(k);
end

Qcond = ks*S*(B(2)-B(1))/(ls_ref/p);

% Inflow/Outflow heat flux

Qi = rho_w*Cp_w*qi*(Ti-Tp) ;

% Rain heat flux

Qr = rho_w*Cp_w*qr*(Ta - Tp)*S;

Tp = Tp + s_t*1/(rho_w*V*Cp_w)*(Qra_p+Qra_s+Qra_a+Qev+Qconv+Qcond+Qi+Qr);

end

