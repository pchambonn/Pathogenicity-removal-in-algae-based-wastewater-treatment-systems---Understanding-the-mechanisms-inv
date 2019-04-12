%% Main pond simulation by step

% This script aims at simulating a pond behaviour  the step by step
% calculation of the different variables.
% This modellling is hoped to offer the advantage of using changeable
% inputs to test for hypothesis of optimization of HRAP for disinfection
% purposes.
% This script is in addition aimed at simulating different optimization
% strategies for HRAP disinfection purpose (see Chambonniere, 2019 PhD
% thesis manuscript.

%% THIS MAIN WAS USED FOR THE VALIDATION OF THE MODELLING RESULTS %%

clearvars

p_t = 100;                                                                                  % Refinement factor on the step of time of weather data inputs 
f_refine = 3;                                                                               % Refinement factor of the calcualtion step of time to be used when TC or DO become negative
t_cut = 1;                                                                                    % step of time under which the computing is cut (s)

%% INPUTS

str = 'Data_Palmerston North';                                                % string of the file name of input weather data 
Time_start = datenum(2016,07,22,10,0,0);                           % date of the first data point (to be entered manually)
str = strcat(str,'_model');

%% IMPORTANT MODELLING PARAMETERS

% l_x: pond length (m)
% l_y: pond width (m2)
% d_obj: design depth (m)
% sedim: sedimentation factor (-)
% PE: algae photosynthetic efficiency (-)
% NO3: HRAP nitrate concentration (mg/L)
% HRT_1: design HRAP HRT (d)
% Kla_O2_1: HRAP aeration rate coefficient (s-1)

l_x = 1;                                                                            
l_y = 3.42;
d_obj = 0.25;
sedim = 0.8;
PE = 0.02;
NO3 = 0;
HRT_1 = 7.9;                                                                         
Kla_O2_1 = 1.0*10^(-4);




%% Hypothese tested
% Hyp: if Hyp == 1: testing normal modelling;
% Hyp = 6: pond operated with no volume control but a constant volume is
    % withdrawn at a constant hour based on HRT. This hour will be based on
    % set_time_1.
% n_hrap: number of HRAP in series tested

n_hrap = 1;                                                                                   
Hyp = 1;
set_time_1 = 18;


%% Simulation Initialisation

Tp_1 = 7.0 + 273.15;                                                                % Initial HRAP broth temperature (K)
X_algae_1 = 0.05;                                                                    % Initial algae concentration (kg/m3)
IC_1 = 22;                                                                                  % Initial inorganic carbon concentration (g C/m3)
DO_1 = 7.8;                                                                               % Initial DO concentration (g O2/m3)
bCOD_1 = 1.5;                                                                          % Initial bCOD concentration (g O2/m3)
X_bacteria_1 = 85;                                                                   % Initial heterotrophic biomass concentration (g/m3)
IP_1 = 0.8;                                                                                 % Initial inorganic phosphorous concentration (g P/m3)
IN_1 = 4;                                                                                     % Initial ammonical nitrogen concentration (g N/m3)
Sigma_1 = 1.4;                                                                         % Initial inert charge balance in the broth (mol eq/m3)
X_nit_1 = 2.64;                                                                          % Initialnitrifying biomass concentration (g/m3)
C_coli_1 = 1.3*10^9;                                                                % Initial E. coli cell count (MPN/m3)

%% Inlet characteristics

T_inlet = 15 + 273.15;                                                              % Wastewater (WW) temperature (K)
COD_IN = 300;                                                                         % WW COD concentration (g O2/m3)
C_IN = 4.7*10^10;                                                                    % WW E. coli cell count (MPN/m3)
IC_IN = 50;                                                                                 % WW inorganic carbon concentration (gC/m3)
Sigma_IN = 2.4;                                                                        % WW inert ions mass balance (mol eq/m3)
IP_IN = 6.12*31/95;                                                                  % WW iorganic phosphorous (g P/m3)
IN_IN = 28.4;                                                                              % WW ammoniacal N (g N/m3) 
X_bacteria_IN = 0;                                                                    % WW heterotrophic biomass concentration (g VSS/m3)

%% Physics parameters

% rho_w : density of pond water (kg/m3)
% Cp_w : speciffic heat capacity of pond water (J/kg/K)
% eps_w : water emissivity
% sigma_stefan : Stefan-Boltzmann constant (W/m²/K) 
% fa : algal absorption fraction (%)
% eps_a : air emissivity
% Lw : water latent heat (J/kg)
% Nu_a : air kinematic viscosity (m²/s)
% R : ideal gas constant (Pa.m3/mol/K)
% Mw : molecular weight of water (kg/mol)
% lambda_a : air thermal conductivity (W/m/k)
% alpha_a : air thermal diffusivity (m²/s)
% Cp_s: soil thermal capacity (J.kg-1.K-1)
% rho_s: soil volumetric mass (kg.m-3)
% ks : soil thermal conductivity (W/m/K)
% alpha_s: soil thermal diffusivity (m2/s)
% Ts_ref: soil reference temperature (K)


rho_w = 998 ;
Cp_w = 4180 ;
eps_w = 0.97;
sigma_stephan = 5.67*10^(-8);
fa = PE;
eps_a = 0.8 ;
Lw = 2.45 * 10^6 ;
Nu_a =1.5 * 10^(-5);
R = 8.314 ;
Mw = 0.018 ;
lambda_a = 2.6 *10^(-2);
alpha_a =2.2*10^(-5);
Cp_s = 1250; 
rho_s = 1900 ;
ks = 1.7;
alpha_s = ks/(Cp_s*rho_s);
r = 0.03;
Ts_ref = 13.6 + 273.15;

%% Chemistry parameters

% KP2: equilibrium constant HPO4/H2PO4
% KP3: equilibrium constant PO4/HPO4
% KN: equilibrium constant NH3/NH4

KP2 =10^(-7.21);
KP3 = 10^(-12.66);
KN = 10^(-9.25);

%% Biological Parameters
    % Algae

K_C = 0.00432;                                                                         % Affinity algal growth with inorganic carbon (g C/m3)
K_DO_algae = 0.02;                                                                % Affinity algal decay with DO (g O2/m3)
HV = 24.7*10^6;                                                                        % Algae heat value (J/kg VSS)
eps_eff = 1.0;                                                                            % Safety factors accoutning for extra energy required by algae in non optimal growth conditions (cf Bechet et al. 2013);
lambda_algae_light = 0.12/24/3600;                                     % rate of respiration of algae in light conditions(kg VSS decayed/kg VSS produced/s)
lambda_algae_dark = 1.2*0.12/24/3600;                             % rate of respiration of algae in dark conditions(kg VSS decayed/kg VSS produced/s)

Y_DO_a_g = 1.5273;                                                               % Utilization rate of DO by algae during growth (g O2/g VSS)
Y_DO_a_d = -1.5273;                                                              % Utilization rate of DO by algae during decay (g O2/g VSS)
Y_C_a_g = -0.5096;                                                                 % Utilization rate of IC by algae during growth (g C/g VSS)
Y_C_a_d = 0.5096;                                                                  % Utilization rate of IC by algae during decay (g C/g VSS)
Y_N_a_g = -0.04145;                                                               % Utilization rate of AIN by algae during growth (g N/g VSS)
Y_N_a_d = 0.04145;                                                                % Utilization rate of AIN by algae during decay (g N/g VSS)

    % Heterotrophic bacteria
    
mu_20 = 6/24/3600;                                                                 % maximal specific growth at 20°C (g VSS/g VSS/s)
Ks_20 = 20;                                                                               % affinity on bCOD at 20°C (mg bCOD/L)
Y_20 = 0.4479;                                                                          % Heterotrophic biomass yield from bCOD consumption (g VSS/g BCOD)
kd_20 = 0.12/24/3600;                                                             %  decay rate at 20 °C (g VSS/g VSS/s)
K_DO_bacteria = 0.2;                                                              % affinity of heterotrophic bacteria on DO (g O2/m3)
teta_mu = 1.07;                                                                         % correction coefficient of max specific growth for temperature
teta_kd = 1.04;                                                                          % correction coefficient of decay rate for temperature
teta_Ks = 1;                                                                               % correction coefficient of affinity on bCOD for temperature
fs_i = 0.15;                                                                                 % debris fraction production (g / g VSS)

Y_DO_b_g = -1.084;                                                                % Utilization rate of DO by heterotrophic bacteria during growth (g O2/g VSS)
Y_C_b_g = 0.2190;                                                                  % Utilization rate of IC by heterotrophic bacteria during growth (g C/g VSS)
Y_N_b_g = -0.03639;                                                               % Utilization rate of AIN by heterotrophic bacteria during growth (g N/g VSS)
Y_DO_b_d = -1.4159;                                                              % Utilization rate of DO by heterotrophic bacteria during decay (g O2/g VSS)
Y_C_b_d = 0.5310;                                                                  % Utilization rate of IC by heterotrophic bacteria during decay (g C/g VSS)
Y_N_b_d = 0.1239;                                                                  % Utilization rate of AIN by heterotrophic bacteria during decay (g N/g VSS)
    

    % Nitrifiers

Y_C_n = -0.5310;                                                                      % Utilization rate of IC by nitrifiers during growth (g C/g VSS)
Y_DO_n = -26.9171;                                                                % Utilization rate of DO by nitrifiers during growth (g O2/g VSS)
Y_alk_n = -0.4425;                                                                    % Utilization rate of inert charges by nitrifiers during growth (mol/gVSS <=> (meq/L)/(gVSS/m3))
Y_N_n = -6.321;                                                                        % Utilization rate of AIN by nitrifiers during growth (g N/g VSS)

    % Miscellaneous
    
sigma_a = 64;                                                                           % transmitance of algal broth (m-1)
alpha_dis = 0.0678;                                                                 % proportionality factor for disinfection rate against sunlight intensity (d-1.m2.J-1)


%% Creation of outputs vectors

Zero = cell(1,n_hrap);

X_algae = Zero;
bCOD = Zero;
X_bacteria = Zero;
C_coli = Zero;
X_debris_a = Zero;
X_debris_b = Zero;
X_nit =  Zero;
DO = Zero; 
pH = Zero; 
H2CO3 = Zero; 
HCO3 = Zero; 
CO3 = Zero; 
IC = Zero;
OH = Zero;
H = Zero;
NH3 = Zero;
NH4 = Zero;
IN = Zero;
H2PO4 = Zero;
HPO4 = Zero;
PO4 = Zero;
IP = Zero;
Sigma = Zero;
T_soil = Zero;
me = Zero;
Tp = Zero;
d = Zero;
HRT = Zero;
Kla_O2 = Zero;
Kla_CO2 = Zero;
Kla_NH3 = Zero;
qi = Zero;
qo = Zero;
kd = Zero;
k_nat = Zero;
k_pH = Zero;
k_sun = Zero;
E_coli_kill_nat = Zero;
E_coli_kill_pH = Zero;
E_coli_kill_sun = Zero;
Time = Zero;

%% Calculations

%% Loop for each hrap of the series

for k_hrap = 1:n_hrap
    
%% Variables initialisation
% Call for weather data
    
    raw_data = xlsread(char(str));
 
% Calculation parameters
    % n : step of time for calculation
    % p : step of depth for soil temperature calculation
    % s_t : step of time (s)

    n = length(raw_data(:,1));
    p = 100;
    s_t = (raw_data(2,1)-raw_data(1,1))*24*3600;

% Weather data are first read and stored as individual vectors
    Hs_dir = raw_data(:,11)';
    Hs_dif = raw_data(:,12)';
    Ta = raw_data(:,3)'+ 273.15;
    Tdp = raw_data(:,4)'+ 273.15;
    RH = zeros(n,1);
    for z =1:n
        RH(z) = Pvap(Tdp(z))/Pvap(Ta(z)) ;
    end
    Nu = raw_data(:,5)';
    qr = raw_data(:,13)/1000/s_t;

% Weather data are then linearly interpolated to obtain data on the initial
% time step divided by the factor p_t
    Hs_dir = linear_interpolation_vector(Hs_dir,p_t);
    Hs_dif = linear_interpolation_vector(Hs_dif,p_t);
    Ta = linear_interpolation_vector(Ta,p_t);
    RH = linear_interpolation_vector(RH,p_t);
    Nu = linear_interpolation_vector(Nu,p_t);
    Dw_a = -2.775*10^(-6)+4.479*10^(-8)*Ta + 1.656*10^(-10)*Ta.^2;
    
    A = zeros(p_t*(n-1) +1,1);
    z = 1;
    for z = 1:n-1
        for k = 1:p_t
            A(p_t*(z-1) + k) = qr(z);
        end
    end
    qr = A;

    n = (n-1)*p_t + 1;
    s_t = s_t/p_t;
    s_t_reference = s_t;
    
    %% Creation of outputs
    
    X_algae{k_hrap} = zeros(n,1);
    bCOD{k_hrap} = zeros(n,1);
    X_bacteria{k_hrap} = zeros(n,1);
    C_coli{k_hrap} = zeros(n,1);
    X_debris_a{k_hrap} = zeros(n,1);
    X_debris_b{k_hrap} = zeros(n,1);
    X_nit{k_hrap} =  zeros(n,1);
    DO{k_hrap} = zeros(n,1); 
    pH{k_hrap} = zeros(n,1); 
    H2CO3{k_hrap} = zeros(n,1); 
    HCO3{k_hrap} = zeros(n,1); 
    CO3{k_hrap} = zeros(n,1); 
    IC{k_hrap} = zeros(n,1);
    OH{k_hrap} = zeros(n,1);
    H{k_hrap} = zeros(n,1);
    NH3{k_hrap} = zeros(n,1);
    NH4{k_hrap} = zeros(n,1);
    IN{k_hrap} = zeros(n,1);
    H2PO4{k_hrap} = zeros(n,1);
    HPO4{k_hrap} = zeros(n,1);
    PO4{k_hrap} = zeros(n,1);
    IP{k_hrap} = zeros(n,1);
    Sigma{k_hrap} = zeros(n,1);
    T_soil{k_hrap} = zeros(n,p);
    me{k_hrap} = zeros(n,1);
    Tp{k_hrap} = zeros(n,1);
    d{k_hrap} = zeros(n,1);
    HRT{k_hrap} = zeros(n,1);
    Kla_O2{k_hrap} = zeros(n,1);
    Kla_CO2{k_hrap} = zeros(n,1);
    Kla_NH3{k_hrap} = zeros(n,1);
    qi{k_hrap} = zeros(n,1);
    qo{k_hrap} = zeros(n,1);
    kd{k_hrap} = zeros(n,1);
    k_nat{k_hrap} = zeros(n,1);
    k_pH{k_hrap} = zeros(n,1);
    k_sun{k_hrap} = zeros(n,1);
    E_coli_kill_nat{k_hrap} = zeros(n,1);
    E_coli_kill_pH{k_hrap} = zeros(n,1);
    E_coli_kill_sun{k_hrap} = zeros(n,1);

    Time{k_hrap} = zeros(n + 1,1);
    
%% Variable initialisation

    Time{k_hrap}(1) = Time_start;
    d{k_hrap}(1) = d_obj;
    HRT{k_hrap}(1) = HRT_1/n_hrap;                                       % underlying hypothesis of equal HRT in all ponds in series
    Kla_O2{k_hrap}(1) = Kla_O2_1;
    Kla_CO2{k_hrap}(1) = Kla_O2{1}(1)*(1.92/2.1)^(1/2);   
    Kla_NH3{k_hrap}(1) = Kla_O2{1}(1)*(1.64/2.1)^(1/2);
    Tp{k_hrap}(1) = Tp_1;
    X_algae{k_hrap}(1) = X_algae_1;
    IC{k_hrap}(1) = IC_1;
    DO{k_hrap}(1) = DO_1;
    bCOD{k_hrap}(1) = bCOD_1;
    X_bacteria{k_hrap}(1) = X_bacteria_1;
    IP{k_hrap}(1) = IP_1;
    IN{k_hrap}(1) = IN_1;
    Sigma{k_hrap}(1) = Sigma_1;
    X_nit{k_hrap}(1) = X_nit_1;
    for k = 1:p
        T_soil{k_hrap}(1,k) = Tp{k_hrap}(1) + (k - 1)/(p - 1) * (Ts_ref - Tp{k_hrap}(1));
    end
    C_coli{k_hrap}(1) = C_coli_1;

    Time{k_hrap}(2) = Time{k_hrap}(1)+ s_t/24/3600;           % The Time vector is incremented at each calculation since the time step is susceptible to change during calculations: . The first increment is performed in the initialization.
    
%% Loop of calculation

    i = 1;
    while i < n

        s_t = (Time{k_hrap}(i+1) - Time{k_hrap}(i))*24*3600; % Calculation of the time step (s) on which the current calculation will be performed             
        if s_t < 0                                                                             % This loop was added to stop the calculation for both erroneous situation (which would trigger infinite loops)
            break                         
        else
            if s_t < t_cut/f_refine
                break
            end
        end
%% Hypothesis loop

        if Hyp == 1
        % HRT and aeration coefficient are here set equal to the initial value.
        % Inlet and outlet are flowing normally (i.e. esp_i and eps_o = 1)       
            HRT{k_hrap}(i) = HRT{k_hrap}(1);
            Kla_O2{k_hrap}(i) = Kla_O2{k_hrap}(1);
            Kla_CO2{k_hrap}(i) = Kla_CO2{k_hrap}(1);
            Kla_NH3{k_hrap}(i) = Kla_NH3{k_hrap}(1);
            eps_i = 1;                                                                      
            eps_o = 1;
            option_vol_corr = 0;
            qi{k_hrap}(i) = eps_i*l_x*l_y*d_obj/HRT{k_hrap}(i)/24/3600;
        else
            if Hyp == 6
            % HRT and aeration coefficient are here set equal to the initial value.
            % Inlet is flowing normally (i.e. esp_i = 1)
            % Outlet is normal for first ponds of the series (eps_o = 1), but nihil for last (eps_o = 0) unless the time is within set_time immediate
            % neighbouhood, in which case eps_o = 1
                if k_hrap < n_hrap
                    HRT{k_hrap}(i) = HRT{k_hrap}(1);
                    Kla_O2{k_hrap}(i) = Kla_O2{k_hrap}(1);
                    Kla_CO2{k_hrap}(i) = Kla_CO2{k_hrap}(1);
                    Kla_NH3{k_hrap}(i) = Kla_NH3{k_hrap}(1);
                    eps_i = 1;
                    eps_o = 1;
                    option_vol_corr = 0;
                    qi{k_hrap}(i) = eps_i*l_x*l_y*d_obj/HRT{k_hrap}(i)/24/3600;
                else
                    HRT{k_hrap}(i) = HRT{k_hrap}(1);
                    Kla_O2{k_hrap}(i) = Kla_O2{k_hrap}(1);
                    Kla_CO2{k_hrap}(i) = Kla_CO2{k_hrap}(1);
                    Kla_NH3{k_hrap}(i) = Kla_NH3{k_hrap}(1);
                    option_vol_corr = 1;
                    eps_i = 1;
                    qi{k_hrap}(i) = eps_i*l_x*l_y*d_obj/HRT{k_hrap}(i)/24/3600;
                    if Time{k_hrap}(i) - floor(Time{k_hrap}(i)) == set_time_1/24
                        eps_o = 1;
                    else
                        if Time{k_hrap}(i) - floor(Time{k_hrap}(i)) < set_time_1/24 &&  Time{k_hrap}(i+1) - floor(Time{k_hrap}(i+1)) > set_time_1/24
                            eps_o =1;
                        else
                            eps_o = 0;
                        end
                    end
                end
            end
        end

        %% HRAP feed characteristics
        % For the first pond, HRAP feed characteristics are the
        % characteristics of the wastewater
        % For the following HRAPs, characteristics are the characteristics
        % of the previous HRAPs, the inlet flow being the outlet flow of
        % the previous HRAP.
        if k_hrap == 1
             COD_IN_s = COD_IN;
             IC_IN_s = IC_IN;
             IP_IN_s = IP_IN;
             IN_IN_s = IN_IN;
             Sigma_IN_s = Sigma_IN;
             T_inlet_s = T_inlet;
             X_bacteria_IN_s = X_bacteria_IN;
             C_IN_s = C_IN;
        else
            n_t_hrap_1 = 1;
            while Time{k_hrap - 1}(n_t_hrap_1) < Time{k_hrap}(i)  % This loop finds out the time corresponding previous pond calendar in current calcualtion since both time vectors are different.
                n_t_hrap_1 = n_t_hrap_1 + 1;
            end
            COD_IN_s = 1/0.9*bCOD{k_hrap - 1}(n_t_hrap_1);
            IC_IN_s = IC{k_hrap - 1}(n_t_hrap_1);
            IP_IN_s = IP{k_hrap - 1}(n_t_hrap_1);
            IN_IN_s = IN{k_hrap - 1}(n_t_hrap_1);
            Sigma_IN_s = Sigma{k_hrap - 1}(n_t_hrap_1);
            T_inlet_s = Tp{k_hrap - 1}(n_t_hrap_1);
            X_bacteria_IN_s = X_bacteria{k_hrap - 1}(n_t_hrap_1);
            C_IN_s = C_coli{k_hrap - 1}(n_t_hrap_1);
            qi{k_hrap}(i) = qo{k_hrap - 1}(n_t_hrap_1);
        end

        % Application of the function 
        [ X_algae_c, bCOD_c , X_bacteria_c, C_coli_c  , X_debris_a_c , X_debris_b_c, X_nit_c , DO_c, pH_c, H2CO3_c , HCO3_c , CO3_c  , IC_c, H_c , OH_c , NH3_c, NH4_c, IN_c, H2PO4_c , HPO4_c , PO4_c , IP_c , Sigma_c  , kd_c , k_nat_c , k_pH_c , k_sun_c , Tp_c , me_c , T_soil_c ,  d_c , qi_c , qo_c , E_coli_kill_nat_c , E_coli_kill_pH_c , E_coli_kill_sun_c ] = ...
                    pond_simulation_step_3_HRAP_series( Hs_dir(i) , Hs_dif(i) , X_algae{k_hrap}(i) , IC{k_hrap}(i) ,DO{k_hrap}(i) , bCOD{k_hrap}(i) , X_bacteria{k_hrap}(i) , IP{k_hrap}(i) ,  Sigma{k_hrap}(i) , Ta(i) , RH(i) , qr(i) , T_soil{k_hrap}(i,:)', Tp{k_hrap}(i), X_debris_a{k_hrap}(i) , X_debris_b{k_hrap}(i), X_nit{k_hrap}(i), C_coli{k_hrap}(i), IN{k_hrap}(i) ...
                    , Dw_a(i) , Nu(i) , rho_w , Cp_w , eps_w , sigma_stephan , fa , eps_a , Lw , Nu_a , R , Mw , lambda_a , alpha_a , alpha_s , Ts_ref , ks ...
                    , K_C , PE , HV , eps_eff , K_DO_algae ,lambda_algae_light , lambda_algae_dark , sigma_a , fs_i , sedim ,  Y_DO_a_g , Y_DO_a_d , Y_C_a_g , Y_C_a_d , Y_N_a_g, Y_N_a_d  ...
                    , Y_20 , mu_20 , Ks_20 , kd_20 , teta_mu , teta_Ks , teta_kd , K_DO_bacteria , NO3 , Y_DO_b_g ,  Y_DO_b_d , Y_C_b_g , Y_C_b_d , Y_N_b_g , Y_N_b_d , Y_DO_n, Y_C_n, Y_N_n, Y_alk_n ...
                    , Kla_O2{k_hrap}(i) , Kla_CO2{k_hrap}(i) , Kla_NH3{k_hrap}(i),  KP2, KP3 , KN ...
                    , COD_IN_s , IC_IN_s , IP_IN_s , IN_IN_s, Sigma_IN_s , T_inlet_s, X_bacteria_IN_s , C_IN_s ...
                    , sigma_a, alpha_dis ...
                    , s_t , l_x , l_y , d{k_hrap}(i), d_obj , HRT{k_hrap}(i), qi{k_hrap}(i), eps_i , eps_o , option_vol_corr  );

            if i == 1                                                                           % In case the inorganic carbon or DO calculated at the first step of time are negative, a message is sent to warn of changing initial conditions since it is not realistic
                if IC_c < 0 || DO_c < 0                                               
                    warning('/!\ Change initial conditions')
                    break
                end
            end

            if IC_c < 0 || DO_c < 0 || Sigma_c < 0                        % In case the IC or DO calculated from the second step of time are negative: if the time step is superior to the accepted limit t_cut: refinement of the time step is applied and the calulation is restarted from the previous step of time; else a default near zero value is applied


                if (Time{k_hrap}(i) - Time{k_hrap}(i-1))*24*3600 > t_cut

                    i = i-1;

                    s_t = (Time{k_hrap}(i+1) - Time{k_hrap}(i))*24*3600;
                    Zeros = zeros( f_refine - 1 , 1);

                    Aux = zeros( f_refine - 1 , 1);                               % The matrix Aux of refined calendar between Time(i) and Time(i+1) is created and inserted into the vector time
                    for k = 1:f_refine - 1
                        Aux(k) = Time{k_hrap}(i) + k*s_t/f_refine/24/3600;
                    end

                    Time{k_hrap} = [ Time{k_hrap}(1:i) ; Aux ; Time{k_hrap}(i+1:end) ];


                % Weather data are linearly interpolated for the newly created date vector
                    Hs_dir = refinement_linear_interpol(Hs_dir,f_refine,i);
                    Hs_dif = refinement_linear_interpol(Hs_dif,f_refine,i);
                    Ta = refinement_linear_interpol(Ta,f_refine,i);
                    RH = refinement_linear_interpol(RH,f_refine,i);
                    Nu = refinement_linear_interpol(Nu,f_refine,i);
                    Dw_a = -2.775*10^(-6)+4.479*10^(-8)*Ta + 1.656*10^(-10)*Ta.^2;
                    qr = [qr(1:i) ; zeros(f_refine - 1,1) ; qr(i+1:end)];

                % Outputs are extended in size by adding a vector of zeros of the size of the added calendar at their end
                    X_algae{k_hrap} = [X_algae{k_hrap} ; Zeros];
                    bCOD{k_hrap} = [bCOD{k_hrap} ; Zeros];
                    X_bacteria{k_hrap} = [X_bacteria{k_hrap} ; Zeros];
                    C_coli{k_hrap} = [C_coli{k_hrap} ; Zeros];
                    X_debris_a{k_hrap} = [X_debris_a{k_hrap} ; Zeros];
                    X_debris_b{k_hrap} = [X_debris_b{k_hrap} ; Zeros];
                    X_nit{k_hrap} = [X_nit{k_hrap} ; Zeros];
                    DO{k_hrap} = [DO{k_hrap} ; Zeros];
                    pH{k_hrap} = [pH{k_hrap} ; Zeros];
                    H2CO3{k_hrap} = [H2CO3{k_hrap} ; Zeros];
                    HCO3{k_hrap} = [HCO3{k_hrap} ; Zeros];
                    CO3{k_hrap} = [CO3{k_hrap} ; Zeros];            
                    IC{k_hrap} = [IC{k_hrap} ; Zeros];
                    OH{k_hrap} = [OH{k_hrap} ; Zeros];
                    H{k_hrap} = [H{k_hrap} ; Zeros];
                    NH3{k_hrap} = [NH3{k_hrap}; Zeros];
                    NH4{k_hrap} = [NH4{k_hrap}; Zeros];
                    IN{k_hrap} = [IN{k_hrap}; Zeros];
                    H2PO4{k_hrap} = [H2PO4{k_hrap} ; Zeros];
                    HPO4{k_hrap} = [HPO4{k_hrap} ; Zeros];
                    PO4{k_hrap} = [PO4{k_hrap} ; Zeros];
                    IP{k_hrap} = [IP{k_hrap} ; Zeros];            
                    Sigma{k_hrap} = [Sigma{k_hrap} ; Zeros];
                    T_soil{k_hrap} = [ T_soil{k_hrap} ; zeros(f_refine - 1 , p)];
                    me{k_hrap} = [me{k_hrap} ; Zeros];
                    Tp{k_hrap} = [Tp{k_hrap} ; Zeros];
                    d{k_hrap} = [d{k_hrap} ; Zeros];
                    qi{k_hrap} = [qi{k_hrap} ; Zeros];
                    qo{k_hrap} = [qo{k_hrap} ; Zeros];
                    kd{k_hrap} = [kd{k_hrap}; Zeros];
                    k_nat{k_hrap} = [k_nat{k_hrap}; Zeros];
                    k_pH{k_hrap} = [k_pH{k_hrap}; Zeros];
                    k_sun{k_hrap} = [k_sun{k_hrap}; Zeros];
                    E_coli_kill_nat{k_hrap} = [E_coli_kill_nat{k_hrap}; Zeros];
                    E_coli_kill_pH{k_hrap} = [E_coli_kill_pH{k_hrap}; Zeros];
                    E_coli_kill_sun{k_hrap} = [E_coli_kill_sun{k_hrap}; Zeros];

                    n = n + f_refine - 1;                                                % The number of computation (i.e. size of output vectors is corrected accordingly)

                    if Time{k_hrap}(i+2) == 0                                     % In case the time for the output of the next calculation doesn't exist (standard case unless a refinement of time step was applied), this time is created based on the standard step of time
                        Time{k_hrap}(i+2) = Time{k_hrap}(i+1) + s_t_reference/24/3600;
                    end
                else                                                                             % We are here in the case where the the time step is below t_cut and the calculations are forced
                    if IC_c < 0
                        IC_c = 0.001;
                    end
                    if DO_c < 0
                        DO_c = 0.001;
                    end
                    if Sigma_c < 0
                        Sigma_c = 0.001;
                    end
                    
                % The outputs vectors are now completed based on the last
                % calculated values. This means there is a known model error
                % since IC_c || DO_c will be negative. The error will larger if
                % t_cut is larger
                    pH{k_hrap}(i) = pH_c;
                    H{k_hrap}(i) = H_c;
                    OH{k_hrap}(i) = OH_c;
                    H2CO3{k_hrap}(i) = H2CO3_c;
                    HCO3{k_hrap}(i) = HCO3_c;
                    CO3{k_hrap}(i) = CO3_c;
                    H2PO4{k_hrap}(i) = H2PO4_c;
                    HPO4{k_hrap}(i) = HPO4_c;
                    PO4{k_hrap}(i) = PO4_c;
                    NH3{k_hrap}(i) = NH3_c;
                    NH4{k_hrap}(i) = NH4_c;
                    Tp{k_hrap}(i+1) = Tp_c;
                    me{k_hrap}(i) = me_c;
                    T_soil{k_hrap}(i+1,:) = T_soil_c;
                    d{k_hrap}(i+1) = d_c;
                    X_algae{k_hrap}(i+1) = X_algae_c;
                    X_bacteria{k_hrap}(i+1) = X_bacteria_c;
                    bCOD{k_hrap}(i+1) = bCOD_c;
                    X_debris_a{k_hrap}(i+1) = X_debris_a_c;
                    X_debris_b{k_hrap}(i+1) = X_debris_b_c;
                    X_nit{k_hrap}(i+1) = X_nit_c;
                    DO{k_hrap}(i+1) = DO_c;
                    IC{k_hrap}(i+1) = IC_c;
                    IP{k_hrap}(i+1) = IP_c;
                    IN{k_hrap}(i+1) = IN_c;
                    Sigma{k_hrap}(i+1) = Sigma_c;
                    C_coli{k_hrap}(i+1) = C_coli_c;
                    kd{k_hrap}(i) = kd_c;
                    k_nat{k_hrap}(i) = k_nat_c;
                    k_pH{k_hrap}(i) = k_pH_c;
                    k_sun{k_hrap}(i) = k_sun_c;
                    qi{k_hrap}(i) = qi_c;
                    qo{k_hrap}(i) = qo_c;
                    E_coli_kill_nat{k_hrap}(i) = E_coli_kill_nat_c;
                    E_coli_kill_pH{k_hrap}(i) = E_coli_kill_pH_c;
                    E_coli_kill_sun{k_hrap}(i) = E_coli_kill_sun_c;
                    
                    if Time{k_hrap}(i+2) == 0                                     % In case the time for the output of the next calculation doesn't exist (standard case unless a refinement of time step was applied), this time is created based on the standard step of time
                        Time{k_hrap}(i+2) = Time{k_hrap}(i+1) + s_t_reference/24/3600;
                    end

                    i = i+1;
                end

            else                                                                                 % Case where IC_c and DO_c were computed normally: the calculation is completed
                % The outputs vectors are now completed based on the
                % calculated values.
                    pH{k_hrap}(i) = pH_c;
                    H{k_hrap}(i) = H_c;
                    OH{k_hrap}(i) = OH_c;
                    H2CO3{k_hrap}(i) = H2CO3_c;
                    HCO3{k_hrap}(i) = HCO3_c;
                    CO3{k_hrap}(i) = CO3_c;
                    H2PO4{k_hrap}(i) = H2PO4_c;
                    HPO4{k_hrap}(i) = HPO4_c;
                    PO4{k_hrap}(i) = PO4_c;
                    NH3{k_hrap}(i) = NH3_c;
                    NH4{k_hrap}(i) = NH4_c;
                    Tp{k_hrap}(i+1) = Tp_c;
                    me{k_hrap}(i) = me_c;
                    T_soil{k_hrap}(i+1,:) = T_soil_c;
                    d{k_hrap}(i+1) = d_c;
                    X_algae{k_hrap}(i+1) = X_algae_c;
                    X_bacteria{k_hrap}(i+1) = X_bacteria_c;
                    bCOD{k_hrap}(i+1) = bCOD_c;
                    X_debris_a{k_hrap}(i+1) = X_debris_a_c;
                    X_debris_b{k_hrap}(i+1) = X_debris_b_c;
                    X_nit{k_hrap}(i+1) = X_nit_c;
                    DO{k_hrap}(i+1) = DO_c;
                    IC{k_hrap}(i+1) = IC_c;
                    IP{k_hrap}(i+1) = IP_c;
                    IN{k_hrap}(i+1) = IN_c;
                    Sigma{k_hrap}(i+1) = Sigma_c;
                    C_coli{k_hrap}(i+1) = C_coli_c;
                    kd{k_hrap}(i) = kd_c;
                    k_nat{k_hrap}(i) = k_nat_c;
                    k_pH{k_hrap}(i) = k_pH_c;
                    k_sun{k_hrap}(i) = k_sun_c;
                    qi{k_hrap}(i) = qi_c;
                    qo{k_hrap}(i) = qo_c;
                    E_coli_kill_nat{k_hrap}(i) = E_coli_kill_nat_c;
                    E_coli_kill_pH{k_hrap}(i) = E_coli_kill_pH_c;
                    E_coli_kill_sun{k_hrap}(i) = E_coli_kill_sun_c;
                    
                if Time{k_hrap}(i+2) == 0                                         % In case the time for the output of the next calculation doesn't exist (standard case unless a refinement of time step was applied), this time is created based on the standard step of time
                    Time{k_hrap}(i+2) = Time{k_hrap}(i+1) + s_t_reference/24/3600;
                end
            end            
        end

    Time{k_hrap} = Time{k_hrap}(1:n);
end

%% Creation of estimates
% time_compliance: cumulated time spent within compliance of microbial
    % discharge guidelines in Palmerston North (d)
% m3_compliance: cumulated water volume spent within compliance of
% microbial discharge guidelines in Palmerston North (m3)
% m3_treated: total water volume treated (m3)

time_compliance = 0;
m3_compliance = 0;
m3_treated = 0;

N_days = floor(Time{n_hrap}(end) - Time{n_hrap}(1)) + 1;
N_compliance = 0;
C_coli_day = zeros(N_days,1);                                                
N_OUT_day = zeros(N_days,1);
V_OUT_day = zeros(N_days,1);
for i = 1:N_days
    k_low = find(i + Time{n_hrap}(1) - 1 <= Time{n_hrap});
    k_high = find(Time{n_hrap} <= i + Time{n_hrap}(1));
    k = intersect(k_low,k_high);                                                 % vector of the indexes of the day 'i'
    if max(k) == n
        k = k(1:end -1);
    end
    N_OUT_day(i) = sum(qo{n_hrap}(k).*C_coli{n_hrap}(k).*(Time{n_hrap}(k+1) - Time{n_hrap}(k))*24*3600); % Number of E. coli cells reaching the outlet on a given day 'i', calculated from the integration of qo*C_coli (MPN)
    V_OUT_day(i) = sum(qo{n_hrap}(k).*(Time{n_hrap}(k+1) - Time{n_hrap}(k))*24*3600); % Volume of water treated on a given day 'i' (m3)
    if Hyp == 1 
        C_coli_day(i) = mean(C_coli{n_hrap}(k));                     % With continuous operation, the average concentration of E. coli in the outlet is the average of the values over the day
    else
        C_coli_day(i) = N_OUT_day(i)/V_OUT_day(i);             % With non-continuous operation, the average concentration of E. coli in the outlet is the ratio of the E. coli reaching the outlet by the total volume treated over the day
    end
    if C_coli_day(i) < 3*10^8                                                     % the number of days of compliance can now be calculated easily 
        N_compliance = N_compliance + 1;
    end
end

N_IN = 0;
N_OUT = 0;
E_coli_decay_rate_inst_int = zeros(n-1,1);
for i = 1:n - 1
    if C_coli{n_hrap}(i) <= 3*10^8
       time_compliance = time_compliance + Time{n_hrap}(i+1) - Time{n_hrap}(i);
       m3_compliance = m3_compliance + qo{n_hrap}(i)*(Time{n_hrap}(i+1) - Time{n_hrap}(i))*24*3600;
    end
    m3_treated = m3_treated + qo{n_hrap}(i)*(Time{n_hrap}(i+1) - Time{n_hrap}(i))*24*3600;
    N_IN = N_IN + qi{n_hrap}(i)*C_IN*(Time{n_hrap}(i+1) - Time{n_hrap}(i))*24*3600;
    N_OUT = N_OUT + qo{n_hrap}(i)*C_coli{n_hrap}(i)*(Time{n_hrap}(i+1) - Time{n_hrap}(i))*24*3600;
    E_coli_decay_rate_inst_int(i) = kd{n_hrap}(i)*(Time{n_hrap}(i+1) - Time{n_hrap}(i));
end

E_coli_decay_rate_av = 1/(Time{n_hrap}(end) - Time{n_hrap}(1))*sum(E_coli_decay_rate_inst_int); % average E. coli decay rate over the simulation (d-1)

LR_av = log10(N_IN/N_OUT);                                                 % total E. coli log removal in the HRAP over the full simulation
time_compliance;
m3_compliance;
m3_treated;

aa_output = [E_coli_decay_rate_av , LR_av , time_compliance , m3_compliance , m3_treated, m3_compliance/m3_treated*100 , N_compliance];


%% ab_output

Tp_int = zeros(n-1,1);
for i = 1:n - 1
    Tp_int(i) = Tp{n_hrap}(i)*(Time{n_hrap}(i+1) - Time{n_hrap}(i));
end
Tp_av = 1/(Time{n_hrap}(end) - Time{n_hrap}(1))*sum(Tp_int); % average pond temperature over teh simulation (K)

pH_int = zeros(n-1,1);
for i = 1:n - 1
    pH_int(i) = pH{n_hrap}(i)*(Time{n_hrap}(i+1) - Time{n_hrap}(i)); % average pond pH over the simulation (-)
end
pH_av = 1/(Time{n_hrap}(end) - Time{n_hrap}(1))*sum(pH_int);

IC_int = zeros(n-1,1);
for i = 1:n - 1
    IC_int(i) = IC{n_hrap}(i)*(Time{n_hrap}(i+1) - Time{n_hrap}(i)); % average pond IC  over the simulation (g C/m3)
end
IC_av = 1/(Time{n_hrap}(end) - Time{n_hrap}(1))*sum(IC_int);

k_nat_int = zeros(n-1,1);
k_pH_int = zeros(n-1,1);
k_sun_int = zeros(n-1,1);
for i = 1:n - 1
    k_nat_int(i) = k_nat{n_hrap}(i)*(Time{n_hrap}(i+1) - Time{n_hrap}(i));
    k_pH_int(i) = k_pH{n_hrap}(i)*(Time{n_hrap}(i+1) - Time{n_hrap}(i));
    k_sun_int(i) = k_sun{n_hrap}(i)*(Time{n_hrap}(i+1) - Time{n_hrap}(i));
end
k_nat_av = 1/(Time{n_hrap}(end) - Time{n_hrap}(1))*sum(k_nat_int);  % average decay rate due to natural die-off over the simulation
k_pH_av = 1/(Time{n_hrap}(end) - Time{n_hrap}(1))*sum(k_pH_int);   % average decay rate due to pH toxicity over the simulation
k_sun_av = 1/(Time{n_hrap}(end) - Time{n_hrap}(1))*sum(k_sun_int); % average decay rate due to sunlight mediated disinfection over the simulation

ab_output = [pH_av , Tp_av - 273.15 , IC_av , k_nat_av , k_pH_av , k_sun_av];

%% ab_output summer
% the same out as previously described in ab_output are isolated over summer

k_summer = find(Time{n_hrap} >= datenum(2016,12,01) & Time{n_hrap} < datenum(2017,03,01));
n_summer = length(k_summer);
    
Tp_int_summer = zeros(n_summer-1,1);
for i = 1:n_summer - 1
    Tp_int_summer(i) = Tp{n_hrap}(k_summer(i))*(Time{n_hrap}(k_summer(i+1)) - Time{n_hrap}(k_summer(i)));
end
Tp_av_summer = 1/(Time{n_hrap}(k_summer(end)) - Time{n_hrap}(k_summer(1)))*sum(Tp_int_summer);

pH_int_summer = zeros(n_summer-1,1);
for i = 1:n_summer - 1
    pH_int_summer(i) = pH{n_hrap}(k_summer(i))*(Time{n_hrap}(k_summer(i+1)) - Time{n_hrap}(k_summer(i)));
end
pH_av_summer = 1/(Time{n_hrap}(k_summer(end)) - Time{n_hrap}(k_summer(1)))*sum(pH_int_summer);


k_nat_int_summer = zeros(n_summer-1,1);
k_pH_int_summer = zeros(n_summer-1,1);
k_sun_int_summer = zeros(n_summer-1,1);
for i = 1:n_summer - 1
    k_nat_int_summer(i) = k_nat{n_hrap}(k_summer(i))*(Time{n_hrap}(k_summer(i+1)) - Time{n_hrap}(k_summer(i)));
    k_pH_int_summer(i) = k_pH{n_hrap}(k_summer(i))*(Time{n_hrap}(k_summer(i+1)) - Time{n_hrap}(k_summer(i)));
    k_sun_int_summer(i) = k_sun{n_hrap}(k_summer(i))*(Time{n_hrap}(k_summer(i+1)) - Time{n_hrap}(k_summer(i)));
end
k_nat_av_summer = 1/(Time{n_hrap}(k_summer(end)) - Time{n_hrap}(k_summer(1)))*sum(k_nat_int_summer);
k_pH_av_summer = 1/(Time{n_hrap}(k_summer(end)) - Time{n_hrap}(k_summer(1)))*sum(k_pH_int_summer);
k_sun_av_summer = 1/(Time{n_hrap}(k_summer(end)) - Time{n_hrap}(k_summer(1)))*sum(k_sun_int_summer);

ab_output_summer = [pH_av_summer , Tp_av_summer - 273.15 , k_nat_av_summer , k_pH_av_summer , k_sun_av_summer];

%% ab_output_winter
% the same out as previously described are isolated over winter

k_winter_1 = find(Time{n_hrap} >= datenum(2016,08,01) & Time{n_hrap} < datenum(2016,09,01));
n_winter_1 = length(k_winter_1);
    
Tp_int_winter_1 = zeros(n_winter_1-1,1);
for i = 1:n_winter_1 - 1
    Tp_int_winter_1(i) = Tp{n_hrap}(k_winter_1(i))*(Time{n_hrap}(k_winter_1(i+1)) - Time{n_hrap}(k_winter_1(i)));
end
Tp_av_winter_1 = 1/((Time{n_hrap}(k_winter_1(end)) - Time{n_hrap}(k_winter_1(1))))*sum(Tp_int_winter_1);

pH_int_winter_1 = zeros(n_winter_1-1,1);
for i = 1:n_winter_1 - 1
    pH_int_winter_1(i) = pH{n_hrap}(k_winter_1(i))*(Time{n_hrap}(k_winter_1(i+1)) - Time{n_hrap}(k_winter_1(i)));
end
pH_av_winter_1 = 1/((Time{n_hrap}(k_winter_1(end)) - Time{n_hrap}(k_winter_1(1))))*sum(pH_int_winter_1);


k_nat_int_winter_1 = zeros(n_winter_1-1,1);
k_pH_int_winter_1 = zeros(n_winter_1-1,1);
k_sun_int_winter_1 = zeros(n_winter_1-1,1);
for i = 1:n_winter_1 - 1
    k_nat_int_winter_1(i) = k_nat{n_hrap}(k_winter_1(i))*(Time{n_hrap}(k_winter_1(i+1)) - Time{n_hrap}(k_winter_1(i)));
    k_pH_int_winter_1(i) = k_pH{n_hrap}(k_winter_1(i))*(Time{n_hrap}(k_winter_1(i+1)) - Time{n_hrap}(k_winter_1(i)));
    k_sun_int_winter_1(i) = k_sun{n_hrap}(k_winter_1(i))*(Time{n_hrap}(k_winter_1(i+1)) - Time{n_hrap}(k_winter_1(i)));
end
k_nat_av_winter_1 = 1/((Time{n_hrap}(k_winter_1(end)) - Time{n_hrap}(k_winter_1(1))))*sum(k_nat_int_winter_1);
k_pH_av_winter_1 = 1/((Time{n_hrap}(k_winter_1(end)) - Time{n_hrap}(k_winter_1(1))))*sum(k_pH_int_winter_1);
k_sun_av_winter_1 = 1/((Time{n_hrap}(k_winter_1(end)) - Time{n_hrap}(k_winter_1(1))))*sum(k_sun_int_winter_1);

k_winter_2 = find(Time{n_hrap} >= datenum(2017,06,01) & Time{n_hrap} < datenum(2017,09,01));
n_winter_2 = length(k_winter_2);
    
Tp_int_winter_2 = zeros(n_winter_2-1,1);
for i = 1:n_winter_2 - 1
    Tp_int_winter_2(i) = Tp{n_hrap}(k_winter_2(i))*(Time{n_hrap}(k_winter_2(i+1)) - Time{n_hrap}(k_winter_2(i)));
end
Tp_av_winter_2 = 1/(((Time{n_hrap}(k_winter_2(end)) - Time{n_hrap}(k_winter_2(1)))))*sum(Tp_int_winter_2);

pH_int_winter_2 = zeros(n_winter_2-1,1);
for i = 1:n_winter_2 - 1
    pH_int_winter_2(i) = pH{n_hrap}(k_winter_2(i))*(Time{n_hrap}(k_winter_2(i+1)) - Time{n_hrap}(k_winter_2(i)));
end
pH_av_winter_2 = 1/(((Time{n_hrap}(k_winter_2(end)) - Time{n_hrap}(k_winter_2(1)))))*sum(pH_int_winter_2);


k_nat_int_winter_2 = zeros(n_winter_2-1,1);
k_pH_int_winter_2 = zeros(n_winter_2-1,1);
k_sun_int_winter_2 = zeros(n_winter_2-1,1);
for i = 1:n_winter_2 - 1
    k_nat_int_winter_2(i) = k_nat{n_hrap}(k_winter_2(i))*(Time{n_hrap}(k_winter_2(i+1)) - Time{n_hrap}(k_winter_2(i)));
    k_pH_int_winter_2(i) = k_pH{n_hrap}(k_winter_2(i))*(Time{n_hrap}(k_winter_2(i+1)) - Time{n_hrap}(k_winter_2(i)));
    k_sun_int_winter_2(i) = k_sun{n_hrap}(k_winter_2(i))*(Time{n_hrap}(k_winter_2(i+1)) - Time{n_hrap}(k_winter_2(i)));
end
k_nat_av_winter_2 = 1/(((Time{n_hrap}(k_winter_2(end)) - Time{n_hrap}(k_winter_2(1)))))*sum(k_nat_int_winter_2);
k_pH_av_winter_2 = 1/(((Time{n_hrap}(k_winter_2(end)) - Time{n_hrap}(k_winter_2(1)))))*sum(k_pH_int_winter_2);
k_sun_av_winter_2 = 1/(((Time{n_hrap}(k_winter_2(end)) - Time{n_hrap}(k_winter_2(1)))))*sum(k_sun_int_winter_2);


ab_output_winter = [(pH_av_winter_1+pH_av_winter_2)/2 , (Tp_av_winter_1+Tp_av_winter_2)/2 - 273.15 , (k_nat_av_winter_1+k_nat_av_winter_2)/2 , (k_pH_av_winter_1+k_pH_av_winter_2)/2 , (k_sun_av_winter_1+k_sun_av_winter_2)/2];

%% ac_output
% this output gather the sum of all E. coli cells removed through naturatl
% disinfection, pH toxicity, or sunlight mediated disinfection

ac_output = zeros(n_hrap,3);
for i = 1:n_hrap
    ac_output(i,:) = [sum(E_coli_kill_nat{i}),sum(E_coli_kill_pH{i}),sum(E_coli_kill_sun{i})];
end



%% Figures
for k = 1:n_hrap
    for j = 1:length(Time{k})
        if (Time{k}(j)>= datenum(2016,12,12) && (Time{k}(j) <= datenum(2016,12,16)) || (Time{k}(j) >= datenum(2016,12,27)) && Time{k}(j) <= datenum(2016,12,31))
            Tp{k}(j) = NaN;
            pH{k}(j) = NaN;
        end
    end
end

figure(100), clf, hold on
lgd = [];
for k = 1:n_hrap
    plot(Time{k},Tp{k}-273.15,'LineWidth',1,'Color',[1-1/(0.5*k + 0.5) 1-1/(0.5*k + 0.5) 1-1/(0.5*k + 0.5)])
    lgd = [lgd ; strcat('HRAP #',num2str(k))];
end
ylabel({'Pond Temperature'  '(^oC)'},'Fontsize',18);
legend(lgd)
ax = gca; xlim(gca,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14; 
ax.XTickLabelRotation = 45;
title('n = 1','FontSize',22)
% datetick('x','dd-mmm-yyyy'); xlim(gca,[Time(850000) Time(end)]);

figure(101), clf, hold on
for k = 1:n_hrap
    plot(Time{k},pH{k},'LineWidth',1,'Color',[1-1/(0.5*k + 0.5) 1-1/(0.5*k + 0.5) 1-1/(0.5*k + 0.5)])
end
legend(lgd)
ylabel('pH','Fontsize',18); 
ax = gca; datetick('x','dd-mmm-yyyy');
ax.FontSize = 14; ylim([7.0 12])
ax.XTickLabelRotation = 45;
title('n = 1','FontSize',22)


figure(102), clf, hold on
for k = 1:n_hrap
    plot(Time{k},DO{k},'LineWidth',1)
end
ylabel({'Dissolved Oxygen'  '(mg.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(103), clf, hold on
for k = 1:n_hrap
    plot(Time{k},X_algae{k}*1000)
end
ylabel({'Algae concentration' '(g TSS.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(104), clf, hold on
for k = 1:n_hrap
    plot(Time{k},IC{k})
end
ylabel({'Inorganic Carbon' '(mg.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(105), clf, hold on
for k = 1:n_hrap
    plot(Time{k},Sigma{k}*50)
end
ylabel({'Alkalinity' '(mg CaCO_3.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(106), clf, hold on
for k = 1:n_hrap
    plot(Time{k},IP{k})
end
ylabel({'PO_4^3^-' '(mg.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(107), clf, hold on
for k = 1:n_hrap
    plot(Time{k},X_bacteria{k})
end
ylabel({'Heterotrophic bacteria' '(mg.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(108), clf, hold on
for k = 1:n_hrap
    plot(Time{k},bCOD{k})
end
ylabel({'bCOD' '(mg.L^-^1)'},'Fontsize',18);
xlim(gca,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax = gca; xlim(gca,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;


figure(109), clf, hold on
for k = 1:n_hrap
    plot(Time{k},X_nit{k})
end
ylabel({'Nitrifying bacteria' '(mg.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(110), clf, hold on
for k = 1:n_hrap
    plot(Time{k}, X_algae{k}*1000 + (X_nit{k} + X_bacteria{k} + X_debris_a{k})/0.9)
end
ylabel({'TSS concentration' '(mg.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(111), clf, hold on
for k = 1:n_hrap
    plot(Time{k}, IN{k})
end
ylabel({'Ammoniac N' '(mg N.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(112), clf, hold on
for k = 1:n_hrap
    plot(Time{k}, d{k})
end
ylabel({'HRAP depth' '(m)'},'Fontsize',18);
ax = gca; xlim(gca,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(113), clf, hold on
for k = 1:n_hrap
    plot(Time{k},C_coli{k})
end
ylabel({'E. coli' '(MPN·m^-^3)'},'Fontsize',18);
A = Time{n_hrap}./Time{n_hrap}*3*10^8;
plot(Time{n_hrap},A,'--r','LineWidth',2)
ax = gca; 
set(ax,'Yscale','log')
xlim(ax,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(114), clf, hold on
for k = 1:n_hrap
    plot(Time{k},qo{k})
end
ylabel({'Flow rate outlet' '(m^3)'},'Fontsize',18);
ax = gca; 
xlim(ax,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

C_coli_out = zeros(n-1,1);
for i = 1:n
    if qo{n_hrap}(i) == 0
        C_coli_out(i) = NaN;
    else
        C_coli_out(i) = C_coli{n_hrap}(i)/10^4;
    end
end
figure(115), clf, hold on
plot(Time{n_hrap},C_coli_out,'+k','MarkerSize',10)
ylabel({'{\itE. coli}' '(MPN·100mL^-^1)'});
A = Time{n_hrap}./Time{n_hrap}*3*10^4;
plot(Time{n_hrap},A,'--r','LineWidth',2)
ax = gca; 
set(ax,'Yscale','log')
xlim(ax,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 18;
ax.XTickLabelRotation = 45;


Time_day = zeros(N_days,1);
for i = 1:N_days
    Time_day(i) = Time_start + (i - 1);
end
figure(116), clf, hold on
plot(Time_day,C_coli_day/10^4,'+k','MarkerSize',10)
ylabel({'{\itE. coli}' '(MPN·100mL^-^1)'});
A = Time_day./Time_day*3*10^4;
plot(Time_day,A,'--r','LineWidth',2)
ax = gca; 
set(ax,'Yscale','log')
xlim(ax,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 18;
ax.XTickLabelRotation = 45;
if Hyp == 1
    hypothesis = '; Continuous';
else
    hypothesis = strcat('; Semi - Continuous, ', num2str(set_time_1),'h');
end
str = strcat('Palmerston North', '; n=',num2str(n_hrap),'; HRT=',num2str(HRT_1),'d',hypothesis);
title(str)

