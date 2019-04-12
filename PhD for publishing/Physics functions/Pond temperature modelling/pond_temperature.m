function [ Time, Hs_dir,Hs_dif, Tp, s_t, me, qr, Qra_p, Qra_s, Qra_a, Qev, Qcond, Qconv, Qi, Qr ] = pond_temperature( Tp_i, Ti, rho_w, Cp_w, V, eps_w, sigma, S, fa, eps_a, Lw, L, Nu_a, R, Mw, lambda_a, alpha_a,Cp_s, rho_s, ks, qi, str, p_t) 

% This function aims to model the temperature profile of any opaque water
% body of uniform temperature profile according to Bechet et al. (2011)


%% INPUTS

% Tp_i: intial temperature (°C)
% Ti : temperature of inlet (°C)
% rho_w : density of pond water (kg/m3)
% Cp_w : speciffic heat capacity of pond water (J/kg/K)
% V : pond volume (m3)
% eps_w : water emissivity
% sigma : Stephan-Boltzmann constant (W/m²/K) 
% S : pond surface (m²)
% fa : algal absorption fraction (%)
% eps_a : air emissivity
% Lw : water latent heat (J/kg)
% L : characteristic pond length (m)
% Nu_a : air kinematic viscosity (m²/s)
% R : ideal gas constant (Pa.m3/mol/K)
% Mw : molecular weight of water (kg/mol)
% lambda_a : air thermal conductivity (W/m/k)
% alpha_a : air thermal diffusivity (m²/s)
% Cp_s : soil specific heat capacity ((J/kg/K)
% rho_s : soil density (kg/m3)
% ks : soil thermal conductivity (W/m/K)
% qi : inflow rate (m3/s)
% str : file name of the raw data:
% p_t : refinement factor fot the time step (integer)

%% Import the data

raw_data = xlsread(char(str));

%% Parameters of the calculation

% n : step of time for calculation
% p : step of depth for soil temperature calculation
% s_t : step of time (s)

n = length(raw_data(:,1));
p = 100;
s_t = (raw_data(2,1)-raw_data(1,1))*24*3600;
Time = zeros(1,n);
for i = 2:n
    Time(i)=(i-1)*s_t;    
end

alpha_s = ks/(Cp_s*rho_s);

%% Input functions

% Hs : solar radiation (W/m²)
% Ta : air temperature (K)
% RH : relative humidity of the air above the pond surface
% Nu : wind velocity (m/s) (must be measured at the boundary layer height
% or a corrective factor should be used
% qr : rain waterflow (m3/m²/s)
% Ti : water inflow temperatrue
% Dw_a : mass difusion coefficiet of water vapor in air (m²/s)
% Ts_ref : soil temperature at ls_ref (K)
% Tdp : Temperature of dew point (K)


Hs_dir = raw_data(:,11)';
Hs_dif = raw_data(:,12)';
Ta = raw_data(:,3)'+ 273.15;
Tdp = raw_data(:,4)'+ 273.15;
RH = zeros(1,n);

for i =1:n
    RH(i) = Pvap(Tdp(i))/Pvap(Ta(i)) ;
end

Nu = raw_data(:,5)';
qr = raw_data(:,13)/1000/s_t;
Ti = (Ti + 273.15);
Dw_a = -2.775*10^(-6)+4.479*10^(-8)*Ta + 1.656*10^(-10)*Ta.^2;
Ts_ref = 13.6 + 273.15;

if size(qi) == [1,1]
    qi = qi*ones(1,n);
end

%% Interpolation
% As for the data is generally one hour which is significantly higher than
% the time neededfor significant changes in water temeprature, data was
% linearly interpolate the inputs with a length multiplication by a factor
% p_t (i.e. the time step is divided by p_t) 

Hs_dir = linear_interpolation_vector(Hs_dir,p_t);
Hs_dif = linear_interpolation_vector(Hs_dif,p_t);
Ta = linear_interpolation_vector(Ta,p_t);
RH = linear_interpolation_vector(RH,p_t);
Nu = linear_interpolation_vector(Nu,p_t);
qi = linear_interpolation_vector(qi,p_t);

A = zeros(p_t*(n-1) +1,1);

i = 1;
for i = 1:n
    for k = 1:p_t
        A(p_t*(i-1) + k) = qr(i);
    end
end
qr = A;

Time = linear_interpolation_vector(Time,p_t);
Dw_a = linear_interpolation_vector(Dw_a,p_t);


n = (n-1)*p_t + 1;
s_t = s_t/p_t;

%% Variables

Tp = zeros(1,n);
Qra_p = zeros(1,n);
Qra_s = zeros(1,n);
Qra_a = zeros(1,n);
Qev = zeros(1,n);
Qconv = zeros(1,n);
Qcond = zeros(1,n);
Qi = zeros(1,n);
Qr = zeros(1,n);
me = zeros(1,n);

%% Initialisation

Tp(1) = Tp_i + 273.15;

B = zeros(p,n);                                                                           % B est la matrice de Ts où B(i,j) = Ts(zi,tj)
for k = 1:n
    B(p,k) = Ts_ref ;
end

for k = 1:p
    B(k,1) = Tp(1)+(k-1)/(p-1)*(Ts_ref-Tp(1));
end

%% Calculation

for i = 1:n-1
    
    %% Heat sources & losses
    % Qra_p : heat loss through pond radiation (W)
    % Qra_s : heat gain through solar radiation (W)
    % Qra_a : heat gain through air radiation (W)
    % Qev : heat loss through evaporation (W)
    % me : rate of evaporation (kg/s/m²)
    % Qconv : heat exchange through convection (W)
    % Qcond: heat exchange through conduction (W)
    % Qi:  heat exchange through inlet/outlet (W)
     % Qr: heat exchange through rainfall (W)

    % Pond radiation

    Qra_p(i) = - eps_w*sigma*Tp(i)^4*S ;
    
    % Solar radiation

    Qra_s(i) = eps_w*(1-fa)*(Hs_dir(i)+Hs_dif(i))*S ;

    % Air radiation

    Qra_a(i) = eps_w*eps_a*sigma*Ta(i)^4*S ;

    % Evaporation

    me(i) = Rate_ev(L,Dw_a(i),Nu_a,Nu(i),Tp(i),RH(i),R,Ta(i),Mw);
        
    Qev(i) = - me(i)*Lw*S ;

    % Convection

    hconv = conv_coeff(Nu_a,alpha_a,L,Nu(i),lambda_a); 
    Qconv(i) = hconv*(Ta(i)-Tp(i))*S ;
    
    % Conduction
    
    ls_ref = 4400*sqrt(alpha_s) ;
    
    %Résolution de B pour la colonne i+1
    
    A = alpha_s*s_t/((ls_ref/p)^2);

    B(1,i) = Tp(i);

    for k = 2:p-1 %1 % on doit pouvoir améliorer le nombre d'opération en ne prenant que jusqu'à p-i
        B(k,i+1) = A*(B(k+1,i)-2*B(k,i) + B(k-1,i)) + B(k,i);
    end

    Qcond(i) = ks*S*(B(2,i)-B(1,i))/(ls_ref/p);
    
    % Inflow/Outflow heat flux

    Qi(i) = rho_w*Cp_w*qi(i)*(Ti-Tp(i)) ;
    
    % Rain heat flux
    
    Qr(i) = rho_w*Cp_w*qr(i)*(Ta(i)-Tp(i))*S;
    
    % Computation of Tp
    
    Tp(i+1) = Tp(i) + s_t*1/(rho_w*V*Cp_w)*(Qra_p(i)+Qra_s(i)+Qra_a(i)+Qev(i)+Qconv(i)+Qcond(i)+Qi(i)+Qr(i));

end




end
