%% Main pond simulation by step

% This script aims at simulating a pond behaviour using the step by step
% calculation of the different variables.
% This modellling is hoped to offer the advantage of using changeable
% inputs to test for hypothesis of optimization of HRAP for disinfection
% purposes.

% This hypoteses assocaited to this script were described in Chambonniere
% (2019) PhD thesis.


clearvars

p_t = 100;                                                                                  % Refinement factor on the step of time of weather data inputs 
f_refine = 3;                                                                               % Refinement factor of the calcualtion step of time to be used when TC or DO become negative
t_cut = 1;                                                                                    % step of time under which the computing is cut (s)


%% INPUTS

str = 'Palmy_validation_year_2';                                             % string of the file name of input weather data 
Time_start = datenum(2016,07,22,10,0,0);                           % date of the first data point (to be entered manually)

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
NO3 = 58;
HRT_1 = 7.9;                                                                         
Kla_O2_1 = 1.0*10^(-4);

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
Y_alk_n = -0.4425;                                                                    % Utilization rate of inert charges by nitrifiers during growth (mol eq/gVSS <=> (meq/L)/(gVSS/m3))
Y_N_n = -6.321;                                                                        % Utilization rate of AIN by nitrifiers during growth (g N/g VSS)

    % Miscellaneous
    
sigma_a = 64;                                                                           % transmitance of algal broth (m-1)
alpha_dis = 0.0678;                                                                 % proportionality factor for disinfection rate against sunlight intensity (d-1.m2.J-1)


%% CALCULATIONS
    %% Variables initialisation

% Call for weather data
str = strcat(str,'_model');                                                           
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
for z = 1:n-1
    for k = 1:p_t
        A(p_t*(z-1) + k) = qr(z);
    end
end
qr = A;

n = (n-1)*p_t + 1;
s_t = s_t/p_t;
s_t_reference = s_t;


%% Creation of outputs vectors

X_algae = zeros(n,1);
bCOD = zeros(n,1);
X_bacteria = zeros(n,1);
C_coli = zeros(n,1);
X_debris_a = zeros(n,1);
X_debris_b = zeros(n,1);
X_nit =  zeros(n,1);
DO = zeros(n,1); 
pH = zeros(n,1); 
H2CO3 = zeros(n,1); 
HCO3 = zeros(n,1); 
CO3 = zeros(n,1); 
IC = zeros(n,1);
OH = zeros(n,1);
H = zeros(n,1);
NH3 = zeros(n,1);
NH4 = zeros(n,1);
IN = zeros(n,1);
H2PO4 = zeros(n,1);
HPO4 = zeros(n,1);
PO4 = zeros(n,1);
IP = zeros(n,1);
Sigma = zeros(n,1);
T_soil = zeros(n,p);
me = zeros(n,1);
Tp = zeros(n,1);
d = zeros(n,1);
HRT = zeros(n,1);
Kla_O2 = zeros(n,1);
Kla_CO2 = zeros(n,1);
Kla_NH3 = zeros(n,1);
qi = zeros(n,1);
qo = zeros(n,1);
kd = zeros(n,1);
Time = zeros(n + 1,1);

%% Variable initialisation

Time(1) = Time_start;
d(1) = d_obj;
HRT(1) = HRT_1;
Kla_O2(1) = Kla_O2_1;
Kla_CO2(1) = Kla_O2(1)*(1.92/2.1)^(1/2);
Kla_NH3(1) = Kla_O2(1)*(1.64/2.1)^(1/2);
Tp(1) = Tp_1;
X_algae(1) = X_algae_1;
IC(1) = IC_1;
DO(1) = DO_1;
bCOD(1) = bCOD_1;
X_bacteria(1) = X_bacteria_1;
IP(1) = IP_1;
IN(1) = IN_1;
Sigma(1) = Sigma_1;
X_nit(1) = X_nit_1;
for k = 1:p
    T_soil(1,k) = Tp(1) + (k - 1)/(p - 1) * (Ts_ref - Tp(1));
end
C_coli(1) = C_coli_1;

Time(2) = Time(1) + s_t/24/3600;                                           % The Time vector is incremented at each calculation since the time step is susceptible to change during calculations: . The first increment is performed in the initialization.

%% Loop of calculation



i = 1;
while i < n

    s_t = (Time(i+1) - Time(i))*24*3600;                                  % Calculation of the time step (s) on which the current calculation will be performed
    if s_t < 0                                                                                 % This loop was added to stop the calculation for both erroneous situation (which would trigger infinite loops)                                                                  
        break
    else
        if s_t < t_cut/f_refine
            break
        end
    end

    % HRT and aeration coefficient are here set equal to the initial value.
    % This step was needed because in a following evolution, the same model was
    % implemented for varying HRT and aeration coefficient.

    HRT(i) = HRT(1);
    Kla_O2(i) = Kla_O2(1);
    Kla_CO2(i) = Kla_CO2(1);
    Kla_NH3(i) = Kla_NH3(1);

    % Application of the function 

    [ X_algae_c, bCOD_c , X_bacteria_c, C_coli_c  , X_debris_a_c , X_debris_b_c, X_nit_c , DO_c, pH_c, H2CO3_c , HCO3_c , CO3_c  , IC_c, H_c , OH_c , NH3_c, NH4_c, IN_c, H2PO4_c , HPO4_c , PO4_c , IP_c , Sigma_c  , kd_c , Tp_c , me_c , T_soil_c ,  d_c , qi_c , qo_c ] = ...
                pond_simulation_step_3( Hs_dir(i) , Hs_dif(i) , X_algae(i) , IC(i) ,DO(i) , bCOD(i) , X_bacteria(i) , IP(i) ,  Sigma(i) , Ta(i) , RH(i) , qr(i) , T_soil(i,:)', Tp(i), X_debris_a(i) , X_debris_b(i), X_nit(i), C_coli(i), IN(i) ...
                , Dw_a(i) , Nu(i) , rho_w , Cp_w , eps_w , sigma_stephan , fa , eps_a , Lw , Nu_a , R , Mw , lambda_a , alpha_a , alpha_s , Ts_ref , ks ...
                , K_C , PE , HV , eps_eff , K_DO_algae ,lambda_algae_light , lambda_algae_dark , sigma_a , fs_i , sedim ,  Y_DO_a_g , Y_DO_a_d , Y_C_a_g , Y_C_a_d , Y_N_a_g, Y_N_a_d  ...
                , Y_20 , mu_20 , Ks_20 , kd_20 , teta_mu , teta_Ks , teta_kd , K_DO_bacteria , NO3 , Y_DO_b_g ,  Y_DO_b_d , Y_C_b_g , Y_C_b_d , Y_N_b_g , Y_N_b_d , Y_DO_n, Y_C_n, Y_N_n, Y_alk_n ...
                , Kla_O2(i) , Kla_CO2(i) , Kla_NH3(i),  KP2, KP3 , KN ...
                , COD_IN , IC_IN , IP_IN , IN_IN, Sigma_IN , T_inlet, X_bacteria_IN , C_IN ...
                , sigma_a, alpha_dis ...
                , s_t , l_x , l_y , d(i), d_obj , HRT(i) );

    if i == 1                                                                                   % In case the inorganic carbon or DO calculated at the first step of time are negative, a message is sent to warn of changing initial conditions since it is not realistic
        if IC_c < 0 || DO_c < 0
            warning('/!\ Change initial conditions')
            break
        end
    end
    
    if IC_c < 0 || DO_c < 0                                                           % In case the IC or DO calculated from the second step of time are negative: if the time step is superior to the accepted limit t_cut: refinement of the time step is applied and the calulation is restarted from the previous step of time; else a default near zero value is applied      
        if (Time(i) - Time(i-1))*24*3600 > t_cut           
            i = i-1;
            s_t = (Time(i+1) - Time(i))*24*3600;
            Zeros = zeros( f_refine - 1 , 1);                                 
            
            Aux = zeros( f_refine - 1 , 1);
            for k = 1:f_refine - 1
                Aux(k) = Time(i) + k*s_t/f_refine/24/3600;             % The matrix Aux of refined calendar between Time(i) and Time(i+1) is created and inserted into the vector time
            end
            Time = [ Time(1:i) ; Aux ; Time(i+1:end) ];


% Weather data are linearly interpolated for the newly created date vector

            Hs_dir = refinement_linear_interpol(Hs_dir,f_refine,i);
            Hs_dif = refinement_linear_interpol(Hs_dif,f_refine,i);
            Ta = refinement_linear_interpol(Ta,f_refine,i);
            RH = refinement_linear_interpol(RH,f_refine,i);
            Nu = refinement_linear_interpol(Nu,f_refine,i);
            Dw_a = -2.775*10^(-6)+4.479*10^(-8)*Ta + 1.656*10^(-10)*Ta.^2;
            qr = [qr(1:i) ; zeros(f_refine - 1,1) ; qr(i+1:end)];
            
% Outputs are extended in size by adding a vector of zeros of the size of the added calendar at their end

            X_algae = [X_algae ; Zeros];
            bCOD = [bCOD ; Zeros];
            X_bacteria = [X_bacteria ; Zeros];
            C_coli = [C_coli ; Zeros];
            X_debris_a = [X_debris_a ; Zeros];
            X_debris_b = [X_debris_b ; Zeros];
            X_nit = [X_nit ; Zeros];
            DO = [DO ; Zeros];
            pH = [pH ; Zeros];
            H2CO3 = [H2CO3 ; Zeros];
            HCO3 = [HCO3 ; Zeros];
            CO3 = [CO3 ; Zeros];            
            IC = [IC ; Zeros];
            OH = [OH ; Zeros];
            H = [H ; Zeros];
            NH3 = [NH3; Zeros];
            NH4 = [NH4; Zeros];
            IN = [IN; Zeros];
            H2PO4 = [H2PO4 ; Zeros];
            HPO4 = [HPO4 ; Zeros];
            PO4 = [PO4 ; Zeros];
            IP = [IP ; Zeros];            
            Sigma = [Sigma ; Zeros];
            T_soil = [ T_soil ; zeros(f_refine - 1 , p)];
            me = [me ; Zeros];
            Tp = [Tp ; Zeros];
            d = [d ; Zeros];
            qi = [qi ; Zeros];
            qo = [qo ; Zeros];
            kd = [kd; Zeros];

            n = n + f_refine - 1;                                                        % The number of computation (i.e. size of output vectors is corrected accordingly)

            if Time(i+2) == 0                                                           % In case the time for the output of the next calculation doesn't exist (standard case unless a refinement of time step was applied), this time is created based on the standard step of time
                Time(i+2) = Time(i+1) + s_t_reference/24/3600;
            end
        else                                                                                    % We are here in the case where the the time step is below t_cut and the calculations are forced
            if IC_c < 0
                IC_c = 0.001;
            end
            if DO_c < 0
                DO_c = 0.001;
            end
            
            % The outputs vectors are now completed based on the last
            % calculated values. This means there is a known model error
            % since IC_c || DO_c will be negative. The error will larger if
            % t_cut is larger. 
            
            pH(i) = pH_c;
            H(i) = H_c;
            OH(i) = OH_c;
            H2CO3(i) = H2CO3_c;
            HCO3(i) = HCO3_c;
            CO3(i) = CO3_c;
            H2PO4(i) = H2PO4_c;
            HPO4(i) = HPO4_c;
            PO4(i) = PO4_c;
            NH3(i) = NH3_c;
            NH4(i) = NH4_c;
            Tp(i+1) = Tp_c;
            me(i) = me_c;
            T_soil(i+1,:) = T_soil_c;
            d(i+1) = d_c;
            X_algae(i+1) = X_algae_c;
            X_bacteria(i+1) = X_bacteria_c;
            bCOD(i+1) = bCOD_c;
            X_debris_a(i+1) = X_debris_a_c;
            X_debris_b(i+1) = X_debris_b_c;
            X_nit(i+1) = X_nit_c;
            DO(i+1) = DO_c;
            IC(i+1) = IC_c;
            IP(i+1) = IP_c;
            IN(i+1) = IN_c;
            Sigma(i+1) = Sigma_c;
            C_coli(i+1) = C_coli_c;
            kd(i) = kd_c;
            qi(i) = qi_c;
            qo(i) = qo_c;

            if Time(i+2) == 0                                                           % In case the time for the output of the next calculation doesn't exist (standard case unless a refinement of time step was applied), this time is created based on the standard step of time
                Time(i+2) = Time(i+1) + s_t_reference/24/3600;
            end

            i = i+1;
        end
    
    else                                                                                         % Case where IC_c and DO_c were computed normally: the calculation is completed
        % The outputs vectors are now completed based on the
        % calculated values.
        pH(i) = pH_c;
        H(i) = H_c;
        OH(i) = OH_c;
        H2CO3(i) = H2CO3_c;
        HCO3(i) = HCO3_c;
        CO3(i) = CO3_c;
        H2PO4(i) = H2PO4_c;
        HPO4(i) = HPO4_c;
        PO4(i) = PO4_c;
        NH3(i) = NH3_c;
        NH4(i) = NH4_c;
        Tp(i+1) = Tp_c;
        me(i) = me_c;
        T_soil(i+1,:) = T_soil_c;
        d(i+1) = d_c;
        X_algae(i+1) = X_algae_c;
        X_bacteria(i+1) = X_bacteria_c;
        bCOD(i+1) = bCOD_c;
        X_debris_a(i+1) = X_debris_a_c;
        X_debris_b(i+1) = X_debris_b_c;
        X_nit(i+1) = X_nit_c;
        DO(i+1) = DO_c;
        IC(i+1) = IC_c;
        IP(i+1) = IP_c;
        IN(i+1) = IN_c;
        Sigma(i+1) = Sigma_c;
        C_coli(i+1) = C_coli_c;
        kd(i) = kd_c;
        qi(i) = qi_c;
        qo(i) = qo_c;

        if Time(i+2) == 0                                                               % In case the time for the output of the next calculation doesn't exist (standard case unless a refinement of time step was applied), this time is created based on the standard step of time
            Time(i+2) = Time(i+1) + s_t_reference/24/3600;
        end

    end
end            


Time = Time(1:n);                                                                      % Due to the calulation at the rank i+2 of time vector, this vector has an extra entry. This line gets rid of it.

%% Creation of estimates

% In this section, output indicative of model performance are created (max
% pH and max temperature for each day of the calculation, time cumulated
% over pH 10, over pH 9.4, and over 20°C.

max_pH = [];
max_Tp = [];

i = 1;
while i <= n - 1 
    j = 0;
    while (Time(i + j) - Time(i)) < 1 && (i+j) < n
        j = j + 1;
    end
    max_pH = [max_pH , max(pH(i:i + j - 1))];
    max_Tp = [max_Tp , max(Tp(i:i + j - 1))];
    i = i + j;
end

t_pH_10 = 0;
t_Tp_20 = 0;
t_pH_9_4 = 0;

for i = 1:n - 1
    if pH(i) > 10
        t_pH_10 = t_pH_10 + Time(i+1) - Time(i);
    end
    if Tp(i) > 293.15
        t_Tp_20 = t_Tp_20 + Time(i+1) - Time(i);
    end
    if pH(i) > 9.4
        t_pH_9_4 = t_pH_9_4 + Time(i+1) - Time(i);
    end
end

aa_output = [t_Tp_20, t_pH_10 , mean(max_Tp) - 273.15 , mean(max_pH)];


%% Figures

figure(100), clf, hold on
plot(Time,Tp-273.15,'k','LineWidth',1), ylabel({'Pond' 'Temperature'  '(^oC)'},'Fontsize',18);
ax = gca; xlim(gca,[Time(1) Time(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14; 
ax.XTickLabelRotation = 45;

figure(101), clf, hold on
plot(Time,pH,'k','LineWidth',1), ylabel('pH','Fontsize',18);
ax = gca; datetick('x','dd-mmm-yyyy');
ax.FontSize = 14; ylim([6.0 12])
ax.XTickLabelRotation = 45;


figure(102), clf, hold on
plot(Time,DO,'k','LineWidth',1), ylabel({'Dissolved Oxygen'  '(mg.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time(1) Time(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;


figure(103), clf, hold on
plot(Time,X_algae*1000,'k'), ylabel({'Algae concentration' '(g TSS.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time(1) Time(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(104), clf, hold on
plot(Time,IC,'k'), ylabel({'Inorganic Carbon' '(mg.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time(1) Time(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(105), clf, hold on
plot(Time,Sigma*50,'k'), ylabel({'Alkalinity' '(mg CaCO_3.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time(1) Time(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(106), clf, hold on
plot(Time,IP,'k'), ylabel({'PO_4^3^-' '(mg.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time(1) Time(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(107), clf, hold on
plot(Time,X_bacteria,'k'), ylabel({'Heterotrophic bacteria' '(mg.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time(1) Time(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(108), clf, hold on
plot(Time,bCOD,'k'), ylabel({'bCOD' '(mg.L^-^1)'},'Fontsize',18);
xlim(gca,[Time(1) Time(end)]); datetick('x','dd-mmm-yyyy');
ax = gca; xlim(gca,[Time(1) Time(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;


figure(109), clf, hold on
plot(Time,X_nit,'k'), ylabel({'Nitrifying bacteria' '(mg.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time(1) Time(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(110), clf, hold on
plot(Time, X_algae*1000 + (X_nit + X_bacteria + X_debris_a)/0.9,'k'), ylabel({'TSS concentration' '(mg.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time(1) Time(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;

figure(111), clf, hold on
plot(Time, IN), ylabel({'Ammoniac N' '(mg N.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time(1) Time(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14;
ax.XTickLabelRotation = 45;


%% Figures over reduced periods
% This allows zooming on specific periods for the plots of temperature, pH,
% and DO by varying n_1 and n_2 with 1 <= n_1 < n_2 <= n;

n_1 = n - 30000;
n_2 = n-1;

figure(120), clf, hold on
plot(Time,Tp - 273.15,'k','LineWidth',1), ylabel('Pond Temperature (°C)','Fontsize',18); plot(Time_meas,Tp_meas,'ok');
xlabel('July 2017','Fontsize',18)
ax = gca;
ax.FontSize = 14; 
ax.XTickLabelRotation = 45;
datetick('x','dd-mmm-yyyy'); xlim(ax,[Time(n_1) Time(n_2)]);

figure(121), clf, hold on
plot(Time,pH,'k','LineWidth',1), ylabel('pH','Fontsize',18); plot(Time_meas,pH_meas,'ok'); 
xlabel('July 2017','Fontsize',18)
ax = gca; datetick('x','dd-mmm-yyyy');
ax.FontSize = 14; 
ax.XTickLabelRotation = 45;
xlim(gca,[Time(n_1) Time(n_2)]);

figure(122), clf, hold on
plot(Time,DO,'k','LineWidth',1), ylabel({'Dissolved Oxygen'  '(mg.L^-^1)'},'Fontsize',18); plot(Time_meas,DO_meas,'ok');
xlabel('July 2017','Fontsize',18)
ax = gca; datetick('x','dd-mmm-yyyy');
ax.FontSize = 14; 
ax.XTickLabelRotation = 45;
xlim(gca,[Time(n_1) Time(n_2)]);

%% Wake up - simulation is done

load gong.mat;
sound(y)

