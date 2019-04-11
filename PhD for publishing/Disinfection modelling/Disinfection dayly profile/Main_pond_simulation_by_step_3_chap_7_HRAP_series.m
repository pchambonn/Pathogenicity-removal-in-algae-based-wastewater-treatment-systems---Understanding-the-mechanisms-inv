%% Main pond simulation by step

% This script aims at simulating a pond behaviour using the by step
% function. This modellling is hoped to offer the advantage of using
% changeable inputs to test for hypothesis of optimization of HRAP for
% disinfection purposes
% This script is referred to as #2 as it used an alternative correction in
% case of negative DO or TC where the calculation just goes one step back
% (instead of back to the start)

%% THIS MAIN WAS USED FOR THE VALIDATION OF THE MODELLING RESULTS %%

clearvars

p_t = 10;
 
f_refine = 3;
t_cut = 1; % step of time under which the computing is cut (s)
% When comparing with real data, base case is p_t = 100, f_refine = 3,
% t_cut = 0.01;

% Data_Palmerston North ; Data_Alexandra ; Data_Whangarei ; Data_Napier ;
% Data_Mount Cook
str = 'Palmy_validation_year_2';
n_hrap = 1 % number of HRAPs in series
HRT_1 = 7.9% Median value over the year
Hyp = 1
d_obj = 0.25;

%% Hypothese tested
% Hyp: if Hyp == 1: testing normal modelling;
% Hyp = 6: pond operated with no volume control but a constant volume is
% withdrawn at a constant hour based on HRT. This hour will be based on
% set_time_1.

set_time_1 = 18;

%% INPUTS

str = strcat(str,'_model');
Time_start = datenum(2016,07,22,0,0,0);

% IMPORTANT MODELLING PARAMETERS (i.e. calculation implication, design
% parameters of HRAP, unceratain variables to set up)

l_x = 1;
l_y = 3.42;

sedim = 0.8; % experimentally determined
PE = 0.02; % calibration
NO3 = 58%58;% mg/L, based on medium value


% Initialisation

Kla_O2_1 = 1.0*10^(-4);

Tp_1 = 7.0 + 273.15;
X_algae_1 = 0.05;
IC_1 = 22; 
DO_1 = 7.8;
bCOD_1 = 1.5;
X_bacteria_1 = 85;
IP_1 = 0.8; 
IN_1 = 4;
Sigma_1 = 1.4;
X_nit_1 = 2.64;
C_coli_1 = 1.3*10^9;

% Inlet


T_inlet = 15 + 273.15;
COD_IN = 300; % M&E (a bit higher than my median value)
C_IN = 4.7*10^10;
IC_IN = 50; % calcualted based on alkalinity from PNCC data
Sigma_IN = 2.4; % calculated based on inert ions balance (mostly K+, Na+, Cl- and total hardness) from PNCC data
IP_IN = 6.12*31/95; % (g P/m3) based on median value
IN_IN = 28.4; % median value
X_bacteria_IN = 0; % hypothesis

% Physics

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


KP2 =10^(-7.21);
KP3 = 10^(-12.66);
KN = 10^(-9.25);


% Biological constants
    % For algae

K_C = 0.00432;
K_DO_algae = 0.02;
HV = 24.7*10^6;
eps_eff = 1.0;
lambda_algae_light = 0.12/24/3600;
lambda_algae_dark = 1.2*0.12/24/3600;

Y_DO_a_g = 1.5273;
Y_DO_a_d = -1.5273;
Y_C_a_g = -0.5096;
Y_C_a_d = 0.5096;
Y_N_a_g = -0.04145;
Y_N_a_d = 0.04145;

    % For heterotrophic bacteria
    
mu_20 = 6/24/3600;
Ks_20 = 20;
Y_20 = 0.4479; 
kd_20 = 0.12/24/3600;
K_DO_bacteria = 0.2;
teta_mu = 1.07;
teta_kd = 1.04;
teta_Ks = 1;
fs_i = 0.15;

Y_DO_b_g = -1.084;
Y_C_b_g = 0.2190;
Y_N_b_g = -0.03639;
Y_DO_b_d = -1.4159;
Y_C_b_d = 0.5310;
Y_N_b_d = 0.1239;
% Nitrifiers

Y_C_n = -0.5310;
Y_DO_n = -26.9171;
Y_alk_n = -0.4425; % mol/gVSS <=> (meq/L)/(gVSS/m3)
Y_N_n = -6.321; 

sigma_a = 64; % transmitance of algal broth (m-1)
alpha_dis = 0.0678; % proportionality factor for disinfection rate against sunlight intensity (d-1.m2.J-1)


%% Inputs variable by time step

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

% Loop for each hrap

for k_hrap = 1:n_hrap
    
    %% Meteo data
    
    raw_data = xlsread(char(str));
    n = length(raw_data(:,1));
    p = 100;
    s_t = (raw_data(2,1)-raw_data(1,1))*24*3600;


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
    
%% Initialisation

    Time{k_hrap}(1) = Time_start;


    d{k_hrap}(1) = d_obj;
    HRT{k_hrap}(1) = HRT_1/n_hrap; % underlying hypothesis of equal HRT in all ponds in series

    Kla_O2{k_hrap}(1) = Kla_O2_1;
    Kla_CO2{k_hrap}(1) = Kla_O2{1}(1)*(1.67/2.01)^(1/2); % source for diffusivity in water is engineeringtoolbox, values at 20°C
    Kla_NH3{k_hrap}(1) = Kla_O2{1}(1)*(1.50/2.01)^(1/2);

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

    Time{k_hrap}(2) = Time{k_hrap}(1)+ s_t/24/3600;

    refinement_index = 0;

    i = 1;
    while i < n

        s_t = (Time{k_hrap}(i+1) - Time{k_hrap}(i))*24*3600;
        if s_t < 0
            break
        else
            if s_t < t_cut/f_refine
                break
            end
        end
%% Hypothesis loop

    if Hyp == 1
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
            while Time{k_hrap - 1}(n_t_hrap_1) < Time{k_hrap}(i)
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
        
        [ X_algae_c, bCOD_c , X_bacteria_c, C_coli_c  , X_debris_a_c , X_debris_b_c, X_nit_c , DO_c, pH_c, H2CO3_c , HCO3_c , CO3_c  , IC_c, H_c , OH_c , NH3_c, NH4_c, IN_c, H2PO4_c , HPO4_c , PO4_c , IP_c , Sigma_c  , kd_c , k_nat_c , k_pH_c , k_sun_c , Tp_c , me_c , T_soil_c ,  d_c , qi_c , qo_c , E_coli_kill_nat_c , E_coli_kill_pH_c , E_coli_kill_sun_c ] = ...
                    pond_simulation_step_3_HRAP_series( Hs_dir(i) , Hs_dif(i) , X_algae{k_hrap}(i) , IC{k_hrap}(i) ,DO{k_hrap}(i) , bCOD{k_hrap}(i) , X_bacteria{k_hrap}(i) , IP{k_hrap}(i) ,  Sigma{k_hrap}(i) , Ta(i) , RH(i) , qr(i) , T_soil{k_hrap}(i,:)', Tp{k_hrap}(i), X_debris_a{k_hrap}(i) , X_debris_b{k_hrap}(i), X_nit{k_hrap}(i), C_coli{k_hrap}(i), IN{k_hrap}(i) ...
                    , Dw_a(i) , Nu(i) , rho_w , Cp_w , eps_w , sigma_stephan , fa , eps_a , Lw , Nu_a , R , Mw , lambda_a , alpha_a , alpha_s , Ts_ref , ks ...
                    , K_C , PE , HV , eps_eff , K_DO_algae ,lambda_algae_light , lambda_algae_dark , sigma_a , fs_i , sedim ,  Y_DO_a_g , Y_DO_a_d , Y_C_a_g , Y_C_a_d , Y_N_a_g, Y_N_a_d  ...
                    , Y_20 , mu_20 , Ks_20 , kd_20 , teta_mu , teta_Ks , teta_kd , K_DO_bacteria , NO3 , Y_DO_b_g ,  Y_DO_b_d , Y_C_b_g , Y_C_b_d , Y_N_b_g , Y_N_b_d , Y_DO_n, Y_C_n, Y_N_n, Y_alk_n ...
                    , Kla_O2{k_hrap}(i) , Kla_CO2{k_hrap}(i) , Kla_NH3{k_hrap}(i),  KP2, KP3 , KN ...
                    , COD_IN_s , IC_IN_s , IP_IN_s , IN_IN_s, Sigma_IN_s , T_inlet_s, X_bacteria_IN_s , C_IN_s ...
                    , sigma_a, alpha_dis ...
                    , s_t , l_x , l_y , d{k_hrap}(i), d_obj , HRT{k_hrap}(i), qi{k_hrap}(i), eps_i , eps_o , option_vol_corr  );

            if i == 1
                if IC_c < 0 || DO_c < 0
                    warning('/!\ Change initial conditions')
                    break
                end
            end

            if IC_c < 0 || DO_c < 0 || Sigma_c < 0


                if (Time{k_hrap}(i) - Time{k_hrap}(i-1))*24*3600 > t_cut

                    i = i-1;
                    refinement_index = refinement_index + 1;

                    s_t = (Time{k_hrap}(i+1) - Time{k_hrap}(i))*24*3600;
                    Zeros = zeros( f_refine - 1 , 1);

                    Aux = zeros( f_refine - 1 , 1);
                    for k = 1:f_refine - 1
                        Aux(k) = Time{k_hrap}(i) + k*s_t/f_refine/24/3600;
                    end

                    Time{k_hrap} = [ Time{k_hrap}(1:i) ; Aux ; Time{k_hrap}(i+1:end) ];



                    Hs_dir = refinement_linear_interpol(Hs_dir,f_refine,i);
                    Hs_dif = refinement_linear_interpol(Hs_dif,f_refine,i);
                    Ta = refinement_linear_interpol(Ta,f_refine,i);
                    RH = refinement_linear_interpol(RH,f_refine,i);
                    Nu = refinement_linear_interpol(Nu,f_refine,i);
                    Dw_a = -2.775*10^(-6)+4.479*10^(-8)*Ta + 1.656*10^(-10)*Ta.^2;
                    qr = [qr(1:i) ; zeros(f_refine - 1,1) ; qr(i+1:end)];


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

                    n = n + f_refine - 1;

                    if Time{k_hrap}(i+2) == 0
                        Time{k_hrap}(i+2) = Time{k_hrap}(i+1) + s_t_reference/24/3600;
                    end
                else
                    if IC_c < 0
                        IC_c = 0.001;
                    end
                    if DO_c < 0
                        DO_c = 0.001;
                    end
                    if Sigma_c < 0
                        Sigma_c = 0.001;
                    end

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

                    if Time{k_hrap}(i+2) == 0
                        Time{k_hrap}(i+2) = Time{k_hrap}(i+1) + s_t_reference/24/3600;
                    end

                    i = i+1;
                end

            else
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
                    
                if Time{k_hrap}(i+2) == 0
                    Time{k_hrap}(i+2) = Time{k_hrap}(i+1) + s_t_reference/24/3600;
                end

                i = i+1;
                if rem(i,10000) == 0
                    i/n
                end
            end            
        end

    Time{k_hrap} = Time{k_hrap}(1:n);
end

%% Creation of estimates


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
    k = intersect(k_low,k_high);
    if max(k) == n
        k = k(1:end -1);
    end
    N_OUT_day(i) = sum(qo{n_hrap}(k).*C_coli{n_hrap}(k).*(Time{n_hrap}(k+1) - Time{n_hrap}(k))*24*3600);
    V_OUT_day(i) = sum(qo{n_hrap}(k).*(Time{n_hrap}(k+1) - Time{n_hrap}(k))*24*3600);
    if Hyp == 1 
        C_coli_day(i) = mean(C_coli{n_hrap}(k));
    else
        C_coli_day(i) = N_OUT_day(i)/V_OUT_day(i);
    end
    if C_coli_day(i) < 3*10^8
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

E_coli_decay_rate_av = 1/(Time{n_hrap}(end) - Time{n_hrap}(1))*sum(E_coli_decay_rate_inst_int);

LR_av = log10(N_IN/N_OUT);
time_compliance;
m3_compliance;
m3_treated;

aa_output = [E_coli_decay_rate_av , LR_av , time_compliance , m3_compliance , m3_treated, m3_compliance/m3_treated*100 , N_compliance];

load gong.mat;
sound(y)

Tp_int = zeros(n-1,1);
for i = 1:n - 1
    Tp_int(i) = Tp{n_hrap}(i)*(Time{n_hrap}(i+1) - Time{n_hrap}(i));
end
Tp_av = 1/(Time{n_hrap}(end) - Time{n_hrap}(1))*sum(Tp_int);

pH_int = zeros(n-1,1);
for i = 1:n - 1
    pH_int(i) = pH{n_hrap}(i)*(Time{n_hrap}(i+1) - Time{n_hrap}(i));
end
pH_av = 1/(Time{n_hrap}(end) - Time{n_hrap}(1))*sum(pH_int);

IC_int = zeros(n-1,1);
for i = 1:n - 1
    IC_int(i) = IC{n_hrap}(i)*(Time{n_hrap}(i+1) - Time{n_hrap}(i));
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
k_nat_av = 1/(Time{n_hrap}(end) - Time{n_hrap}(1))*sum(k_nat_int);
k_pH_av = 1/(Time{n_hrap}(end) - Time{n_hrap}(1))*sum(k_pH_int);
k_sun_av = 1/(Time{n_hrap}(end) - Time{n_hrap}(1))*sum(k_sun_int);

ab_output = [pH_av , Tp_av - 273.15 , IC_av , k_nat_av , k_pH_av , k_sun_av];

%% ab_output summer
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

%
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
% datetick('x','dd-mmm-yyyy'); xlim(gca,[Time(850000) Time(end)]);


figure(102), clf, hold on
for k = 1:n_hrap
    plot(Time{k},DO{k},'LineWidth',1)
end
ylabel({'Dissolved Oxygen'  '(mg.L^-^1)'},'Fontsize',18);
ax = gca; xlim(gca,[Time{1}(1) Time{1}(end)]); datetick('x','dd-mmm-yyyy');
% ylim(gca,[0 18])
ax.FontSize = 14;
ax.XTickLabelRotation = 45;
% datetick('x','dd-mmm-yyyy'); xlim(gca,[Time(850000) Time(end)]);

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

