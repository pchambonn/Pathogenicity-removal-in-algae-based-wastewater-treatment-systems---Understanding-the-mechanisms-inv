function [ X_algae , bCOD , X_bacteria, C_coli , X_debris_a, X_debris_b , X_nit , DO, pH ,H2CO3 , HCO3 , CO3 , IC, H , OH, NH3,  NH4, IN , H2PO4 , HPO4 , PO4 , IP , Sigma  , kd , Tp , me , T_soil ,  d , qi , qo ] = ...
      pond_simulation_step_3( Hs_dir_i , Hs_dif_i , X_algae_i , IC_i ,DO_i , bCOD_i , X_bacteria_i , IP_i ,  Sigma_i , Ta_i , RH_i , qr_i , T_soil_i, Tp_i , X_debris_a_i, X_debris_b_i, X_nit_i, C_coli_i, IN_i ...
    , Dw_a , Nu , rho_w , Cp_w , eps_w , sigma_stephan , fa , eps_a , Lw , Nu_a , R , Mw , lambda_a , alpha_a , alpha_s , Ts_ref , ks ...
    , K_C , PE , HV , eps_eff , K_DO_algae ,lambda_algae_light , lambda_algae_dark , sigma , fs_i , sedim ,  Y_DO_a_g , Y_DO_a_d , Y_C_a_g , Y_C_a_d , Y_N_a_g, Y_N_a_d  ...
    , Y_20 , mu_20 , Ks_20 , kd_20 , teta_mu , teta_Ks , teta_kd , K_DO_bacteria , NO3 , Y_DO_b_g ,  Y_DO_b_d , Y_C_b_g ,  Y_C_b_d , Y_N_b_g, Y_N_b_d , Y_DO_n ,Y_C_n , Y_N_n, Y_alk_n ...
    , Kla_O2 , Kla_CO2 , Kla_NH3, KP2 , KP3 , KN ...
    , COD_IN , IC_IN , IP_IN , IN_IN, Alk , T_inlet, X_bacteria_IN , C_IN ...
    , sigma_a, alpha_dis ...
    , s_t , l_x , l_y , d_i, d_obj , HRT_i )

 
%% 06/11/2018: Massive changes
%%

% This function aims at performing a full simulation of an HRAP environmental parameters
% and disinfection performance from design parameters, influent characteristics, and meteorological
% data. The calculation starts at a defined situation and gives the value
% for the next given step of time.
% This script intent is to provide a function suitable for short terms
% predictions, as well as using dynamic values for inputs such as kla.

%% INPUTS

    % Design parameters

% l_x : width of the pond (m)
% l_y : length of the pond (m)
% d_i : depth of the pond (m) 
% HRT : Hydraulic Retention Time (d)
% s_t: time difference between the input time and the output time (s)

    % Calculation options:

% option_dis : choice for disinfection modelling: 1 = Marais (1974); 2 = Mayo
%   (1995) ; 3 = Craggs (2003) ; 4 = manual adjustment; 5 = PhD
%   disinfection
% option_pH: choice for the calculation of pH: 1 = the pH is not included in the simulation 
%   (technically it is assuming it is constant = 7) and inorganic carbon is tweaked to be in excess 
%   (play on TCin); 2 = the pH is calculated accounting only for carbonate
%   species equilibria; 3 = the pH is calcualted accounting for carbonate and phosphorous species equilibria

    % Physico-chemical constants

% rho_w : density of pond water (kg/m3)
% Cp_w : speciffic heat capacity of pond water (J/kg/K)
% eps_w : water emissivity
% sigma_stephan : Stephan-Boltzmann constant (W/m²/K) 
% fa : algal absorption fraction (%)
% eps_a : air emissivity
% Lw : water latent heat (J/kg)
% Nu_a : air kinematic viscosity (m²/s)
% R : ideal gas constant (Pa.m3/mol/K)
% Mw : molecular weight of water (kg/mol)
% lambda_a : air thermal conductivity (W/m/k)
% alpha_a : air thermal diffusivity (m²/s)
% Dw_a: mass diffusion coefficient of water vapor in air (m2/s)
% alpha_s: soil diffusivity (m2/s)
% Ts_ref: soil reference temperature (K)
% ks : soil thermal conductivity (W/m/K)
% r : reflexivity of the pond (%)
% Kla_O2: mass transfer coef of oxygen from the atm to the pond (s-1)
% Kla_CO2 : mass transfer coef of CO2 form atmosphere to water (s-1)
% Alk : water alkalinity (mol eq/L)
% KP1 : equilibrium constant H3PO4 - H2PO4 (source = wiki)
% KP2 : equilibrium constant H2PO4 - HPO4
% KP3: equilibrium constant HPO4 - PO4
% sedim : sedimentation coefficient accounting for the non ideal mixing (no
% units)

    % Biological constants

    %Algae associated
% K_C : affinity of algae on carbon (on Total Carbon mg-C/L, 0.00432 from Solimenio et al. ; 0.035 from Bai et al.)
% K_DO_algae : affinity of algae on DO during respiration (mg/L) (0.02 from Solimeno et al.)
% sigma : optical density of the algal medium (m2/kg)
% PE: algae photosynthetic efficiency
% HV: Heat value of algae (J/kg DW)
% eps_eff: safety factors accoutning for extra energy required by algae in non optimal growth conditions (cf Bechet 2013);
% lambda_algae_light/dark: rate of respiration of algae (kg TSS decayed/kg TSS produced/s). 1 value in the dark, one value in the light; 
% Y_DO_a_g: Yield of DO production associated with algal growth (g DO/g TSS produced)
% Y_DO_a_d: Yield of DO consumption associated with algal decay (g DO/g TSS decayed)
% Y_IC_a_g: Yield of inorganic carbon consumption associated with algal growth (g C/g TSS decayed)
% Y_IC_a_d: Yield of inorganic carbon release associated with algal decay (g C/g TSS decayed)
% Y_alk_a_g: Yield of alkalinity by algal growth (mol eq/g TSS)
    
    %Bacteria associated
% mu_20 : maximal specific growth at 20*C (g VSS/g VSS/s)
% Ks_20 : affinity on bCOD at 20*C (mg bCOD/L)
% Y_20 : Yield bacterial growth at at 20*C (g VSS/g bCOD)
% kd_20 : decay rate (g VSS/g VSS/s)
% teta_mu : correction coefficient of max specific growth for temperature
% teta_kd : correction coefficient of decay rate for temperature
% teta_Ks : correction coefficient of affinity on bCOD for temperature
% fs_i : production share of inerts from decay (rest is bCOD)
% K_DO_bacteria : affinity on dissolved oxygen(mg/L)
% Y_DO_growth : dissolved oxygen consuption associated to bacterial growth
% (g O2/gVSS)
% Y_C_growth : inorganic carbon consumption associated to bacterial growth
% (g C/gVSS)
% Y_DO_decay : dissolved oxygen consumed associated to bacterial decay (g
% O2/g VSS)
% Y_C_decay : Inorganic carbon production associated to bacterial decay (g
% C-CO2/g VSS)
% Y_alk_growth : alkalinity contribution through bacterial growth (mol eq/g
% VSS)
% fs_i: fraction of debris created as inert
% NO3 : NO3 concentration due to nitrifier activity (mg NO3/L)
% Y_O2_nit: DO consumption from nitrifiers activity (g DO/g VSS)
% Y_C_nit: Inorganic carbon consumption from nitrifiers activity (g C/g VSS)
% Y_alk_nit: Alkalinity consumption from nitrifiers activity (meq /g VSS)


    % Inlet characteristics

% T_inlet : Temeprature of inlet (C)
% COD_IN: influent chemical oxygen demand (mg/L)
% C_IN : concentration of influent in total coliforms supposed constant(MPN/m3)
% TC_IN : Total inorganic carbon in influent (mg C/L)
% TP_IN : Total inorganic phosphorous in influent (mg P/L)
% X_bacteria_IN : Heterotrophic bacteria concentration in the inlet (mg/L)


% Initial values of variables to compute

% Hs_dir_i: Sunlight direct radiation (W/m2) 
% Hs_dif_i: Sunlight indirect radiation (W/m2)
% X_algae_i: algal concentration (kg/m3)
% TC_i: Total inorganic carbon concentration (g-C/m3)
% DO_i: Dissolved oxygen cocnentration (mg/L)
% bCOD_i: biodegradable Chemical oxygen demand (g/m3)
% X_bacteria_i: heterotrophic bacteria concentration (g/m3)
% TP_i: Total inorganic phosphorous (g/m3)
% Gamma_i: pond alkalinity (meq/L)
% Ta_i: air temperature (K)
% RH_i: Relative humidity
% T_soil_i: soil temperature profile (K)
% qr_i: rainflow rate (m3/s)
% Nu_i: wind speed (m/s)
 
    % Manual disinfection rate
    
% k_s_m: solar dependent manually adjusted coefficient (m2.W-1.d-1)
% k_pH: pH dependent manually adjusted coefficient (d-1)
% k_20_m: temperature dependent manually adjusted coefficient (d-1)
% teta_m: temperature dependent Arrheniusish manually adjusted coefficient (-)
% k_DO_m: Disolved oxygen concentration dependent manually adjusted coefficient (L.mgDO-1.d-1)



%% OUTPUTS: need to double check all the units)

% X_algae: algae concentration (kg/m3)
% X_bacteria: heterotroph concentration (g/m3)
% bCOD: biodegradable Chemical oxygen demand (g/m3)
% X_debris: debris fraction of the VSS (g/m3)
% X_nit: Nitrifier concentration (g/m3)
% DO: dissolved oxygen concentration (mg/L)
% TC: Total inorganic carbon (mg/L)
% TP: total inorganic phosphorous (mg/L)
% Gamma: Pond alklinity (meq/L)
% C_coli: E. coli concentration (MPN/m3)
% kd: E. coli decay rate (d-1) 
% pH
% Tp: pond temperature (K)
% me: evaporation rate (kg/m2/s)
% T_soilL: soil temperature profile (K)
% H: H+ concnetration (mg/L)
% OH: OH- concentration (mg/L)
% H2CO3: carbonic acid concentration (mg-C/L)
% HCO3: bicarbonate concentration (mg-C/L)
% CO3: carbonate concentration (mg-C/L)
% H3PO4: (mg/L)
% H2PO4: (mg/L)
% HPO4: (mg/L)
% PO4: (mg/L)
% d: pond new depth (m)




%% Function


% PRIMO CALCULATIONS

% S: pond surface (m2)
% V: pond volume (m3)
% qi: inlet flow rate (m3/s)
% qo: outlet flow rate (m3/s)
% bCOD_IN : influent biological chemical oxygen demand (mg/L)

S = l_x*l_y;
V = S*d_i;
qi = V/HRT_i/24/3600;
bCOD_IN = 0.9*COD_IN;

pH = pH_pond_4( IC_i , IP_i, IN_i, Sigma_i/1000, Tp_i, KP2, KP3, KN );  

% Concentration of species in equilibrium calculations

KC1 = K_carbonate_1(Tp_i);
KC2 = K_carbonate_2(Tp_i);
Ke = K_ionic_product(Tp_i);

nh = 10^(-pH);
H = nh*1000; % conversion for mol-H/L to g-H/m3
OH = Ke/nh*1000;

H2CO3 = IC_i/(1 + KC1/nh + KC1*KC2/nh^2);
HCO3 = H2CO3*KC1/nh;
CO3 = H2CO3*KC1*KC2/(nh^2);

H2PO4 = IP_i/(1 + KP2/nh + KP2*KP3/nh);
HPO4 = H2PO4*KP2/nh;
PO4 = H2PO4*KP2*KP3/(nh^2);

NH4 = IN_i/(1 + KN/nh);
NH3 = NH4*KN/nh;

%% Pond Temperature

[Tp, me, T_soil,~,~,~,~,~,~,~ ] = pond_temperature_instant(Tp_i,Ta_i,Hs_dir_i,Hs_dif_i,RH_i,qr_i,Dw_a,Nu,T_soil_i,T_inlet,qi,s_t,V,S,l_y,rho_w,Cp_w,eps_w,sigma_stephan,fa,eps_a,Nu_a,R,Mw,Lw,lambda_a,alpha_a,alpha_s,Ts_ref,ks);
%  Tp = 1.1*Tp;                 
%% Volume variation
qe = me/rho_w*S;
diff = (qi + qr_i - qe)*s_t;
DF = (d_i + (qr_i - qe)*s_t/S)/d_i;

    if diff > 0
        if d_i < d_obj
            d = min(d_obj , d_i + qi*s_t/S + qr_i*s_t/S - qe*s_t/S);
            qo = max( 0 , qi + qr_i - qe - (d_obj - d_i)*S/s_t);
        else
            d = d_obj;
            qo = qi + qr_i - qe;
        end
    else
        d = d_i + qi*s_t/S + qr_i*s_t/S - qe*s_t/S;
        qo = 0;
    end

    if d < 0.8*d_obj
        V_add = (d_obj - d) * S; % leave the possibility here to take a volume added as an ouput and compute the water footprint of the algal pond. Note the addition criteria can be changed, could be used easily in the input.
        d = d_obj;
    end
    
%% Bilogical variables    
    
    
    
% Algae concentration evolution
[ P_algae_gross, P_algae_resp, ~ ] = algae_simple_productivity( Hs_dir_i, Hs_dif_i, Tp_i , sigma, X_algae_i , d_i, S, IC_i, K_C, PE, HV, eps_eff,DO_i,K_DO_algae,lambda_algae_light,lambda_algae_dark,teta_kd);

P_algae_resp = (1 - fs_i)*P_algae_resp;
Q_algae_debris = fs_i*P_algae_resp;
P_algae_net = P_algae_gross - P_algae_resp - Q_algae_debris;

X_algae = X_algae_i/DF + P_algae_net*s_t/V - sedim*qo*X_algae_i*s_t/V;
X_debris_a = X_debris_a_i/DF + Q_algae_debris*s_t/V - sedim*qo*X_debris_a_i*s_t/V; 

% Bacteria concentration evolution
[ Q_bacteria_gross, Q_bacteria_decay ] = heterotrophic_bacteria_growth( Tp_i,mu_20,Ks_20,kd_20,teta_mu,teta_Ks,teta_kd,DO_i,K_DO_bacteria,bCOD_i,X_bacteria_i,d_i,S );
Q_bacteria_debris = fs_i*Q_bacteria_decay;
Q_bacteria_decay = (1-fs_i)*Q_bacteria_decay;    

X_bacteria = X_bacteria_i/DF + qi*s_t*X_bacteria_IN/V - sedim*qo*s_t*X_bacteria_i/V + Q_bacteria_gross*s_t/V - Q_bacteria_decay*s_t/V - Q_bacteria_debris*s_t/V;
X_debris_b = X_debris_b_i/DF + Q_algae_debris*s_t/V - sedim*qo*X_debris_b_i*s_t/V; 

% bCOD concentration evolution
bCOD = bCOD_i/DF - (1/Y_20)*Q_bacteria_gross*s_t/V + qi*bCOD_IN/V*s_t - qo*bCOD_i*s_t/V ;    

% Nitrifiers
P_NO3 = NO3*qo; % NO3 productivity (g NO3/s)
Q_nit = P_NO3*0.0196/0.98*113/62; % g VSS nit/s
X_nit = X_nit_i/DF + Q_nit*s_t/V - sedim*qo*X_nit_i*s_t/V; %neglecting nitrifiers debris and decay


%% Chemical equilibria

% Calculation of DO concentration (mg/L)
   
DO = DO_i + s_t*(saturation_DO_no_salinity(Tp_i) - DO_i)*Kla_O2 + 1000*Y_DO_a_g*P_algae_gross*s_t/V + 1000*Y_DO_a_d*P_algae_resp*s_t/V - qo*s_t*DO_i/V + ... % Impact of inlet is neglected
     + Y_DO_b_g*Q_bacteria_gross*s_t/V + Y_DO_n*Q_nit*s_t/V + Y_DO_b_d*Q_bacteria_decay*s_t/V;

% Calculations of inorganic IC, IP, and IN to prepare the calcualtion of pH:
H2CO3_d_sat = CO2_solubility(Tp_i);
IC_supply = IC_IN*qi*s_t/V + Kla_CO2*(H2CO3_d_sat - H2CO3)*s_t + Y_C_a_g*P_algae_gross*1000*s_t/V + Y_C_a_d*P_algae_resp*1000*s_t/V - qo*IC_i*s_t/V +...
    Y_C_b_g*Q_bacteria_gross*s_t/V + Y_C_n*Q_nit*s_t/V + Y_C_b_d*Q_bacteria_decay*s_t/V; 

IC = IC_i + IC_supply;

IP_supply = -0.01*(P_algae_gross - P_algae_resp)*s_t*1000/V - 0.01*(Q_bacteria_gross - Q_bacteria_decay + Q_nit)*s_t/V + IP_IN*qi*s_t/V - IP_i*qo*s_t/V; 
IP = max(0.01, IP_i + IP_supply);

IN_supply = IN_IN*qi*s_t/V + Kla_NH3*(0 - NH3)*s_t + Y_N_a_g*P_algae_gross*1000*s_t/V + Y_N_a_d*P_algae_resp*1000*s_t/V - qo*IN_i*s_t/V +... % Ammonia saturation concentration is 0 because there is no ammonia in the atmposphere (corroborated by Solimeno et al.)
    Y_N_b_g*Q_bacteria_gross*s_t/V + Y_N_n*Q_nit*s_t/V + Y_N_b_d*Q_bacteria_decay*s_t/V;
IN = max(0.01, IN_i + IN_supply);


% Calculation of inert charges balance evolution

Sigma = Sigma_i + qi*s_t*Alk/V - qo*s_t*Sigma_i/V + Y_alk_n*Q_nit*s_t/V;


%% E. coli disinfection calculations

G_0 = (Hs_dir_i + Hs_dif_i);
Kd = PhD_disinfection(s_t , Tp_i , pH , G_0 , sigma_a , d_i , alpha_dis)/24/3600; % instant decay rate (s-1)
kd = Kd*24*3600;

C_coli = C_coli_i/DF - Kd*C_coli_i*s_t + (C_IN - sedim*C_coli_i)*qi*s_t/V;



end




