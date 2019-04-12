function [ X_algae , bCOD , X_bacteria, C_coli , X_debris_a, X_debris_b , X_nit , DO, pH ,H2CO3 , HCO3 , CO3 , IC, H , OH, NH3,  NH4, IN , H2PO4 , HPO4 , PO4 , IP , Sigma  , kd , k_nat , k_pH , k_sun , Tp , me , T_soil ,  d , qi , qo , E_coli_kill_nat , E_coli_kill_pH , E_coli_kill_sun ] = ...
      pond_simulation_step_3_HRAP_series( Hs_dir_i , Hs_dif_i , X_algae_i , IC_i ,DO_i , bCOD_i , X_bacteria_i , IP_i ,  Sigma_i , Ta_i , RH_i , qr_i , T_soil_i, Tp_i , X_debris_a_i, X_debris_b_i, X_nit_i, C_coli_i, IN_i ...
    , Dw_a , Nu , rho_w , Cp_w , eps_w , sigma_stephan , fa , eps_a , Lw , Nu_a , R , Mw , lambda_a , alpha_a , alpha_s , Ts_ref , ks ...
    , K_C , PE , HV , eps_eff , K_DO_algae ,lambda_algae_light , lambda_algae_dark , sigma , fs_i , sedim ,  Y_DO_a_g , Y_DO_a_d , Y_C_a_g , Y_C_a_d , Y_N_a_g, Y_N_a_d  ...
    , Y_20 , mu_20 , Ks_20 , kd_20 , teta_mu , teta_Ks , teta_kd , K_DO_bacteria , NO3 , Y_DO_b_g ,  Y_DO_b_d , Y_C_b_g ,  Y_C_b_d , Y_N_b_g, Y_N_b_d , Y_DO_n ,Y_C_n , Y_N_n, Y_alk_n ...
    , Kla_O2 , Kla_CO2 , Kla_NH3, KP2 , KP3 , KN ...
    , COD_IN , IC_IN , IP_IN , IN_IN, Alk , T_inlet, X_bacteria_IN , C_IN ...
    , sigma_a, alpha_dis ...
    , s_t , l_x , l_y , d_i, d_obj , HRT_i , qi_i , eps_i , eps_o , option_vol_corr )

 
% This function aims at performing a full simulation of an HRAP
% environmental parameters (temperature, pH, DO, inorganic nutrients,
% bCOD, micro-organisms, and E. coli disinfection)from design parameters,
% influent characteristics, and meteorological data. 
% The calculation starts at a know initial state and gives the value
% for the next step of time, given listed parameters and the duration of
% next calculation. 
% This script in a evolution of pond_simulation_3 adapted to simutions with
% HRAPs in series

%% INPUTS

  % Initial conditions and weather conditions

% Hs_dir_i: Sunlight direct radiation (W/m2) 
% Hs_dif_i: Sunlight indirect radiation (W/m2)
% X_algae_i: algal concentration (kg/m3)
% IC_i: Total inorganic carbon concentration (g-C/m3)
% DO_i: Dissolved oxygen cocnentration (g O2/m3)
% bCOD_i: biodegradable Chemical oxygen demand (g O2/m3)
% X_bacteria_i: heterotrophic bacteria concentration (g/m3)
% IP_i: Total inorganic phosphorous (g P/m3)
% Sigma_i: pond inert charge balance (meq/L)
% Ta_i: air temperature (K)
% RH_i: Relative humidity
% qr_i: rainflow rate (m3/s)
% T_soil_i: soil temperature profile (K)
% Tp_i: pond temperature (K)
% X_debris_a_i: algae debris concentraiton (kg VSS/m3)
% X_debris_b_i: bacteria debris concentration (g VSS/m3)
% X_nit_i: nitrifier concentration (g VSS/m3)    
% C_coli_i: E. coli concentration
% IN_i: Inorganic ammoniacal N concentration
% Nu: wind speed (m/s)

    % Physical constants

% Dw_a: mass diffusion coefficient of water vapor in air (m2/s)
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
% alpha_s: soil diffusivity (m2/s)
% Ts_ref: soil reference temperature (K)
% ks : soil thermal conductivity (W/m/K)

% Biological constants

    %Algae associated
    
% K_C : affinity of algae on carbon (on Total Carbon mg-C/L, 0.00432 from Solimenio et al. ; 0.035 from Bai et al.)
% PE: algae photosynthetic efficiency
% HV: Heat value of algae (J/kg DW)
% eps_eff: safety factors accoutning for extra energy required by algae in non optimal growth conditions (cf Bechet 2013);
% K_DO_algae : affinity of algae on DO during respiration (mg/L) (0.02 from Solimeno et al.)
% lambda_algae_light/dark: rate of respiration of algae (kg TSS decayed/kg TSS produced/s). 1 value in the dark, one value in the light; 
% sigma : optical density of the algal medium (m2/kg)
% fs_i : production share of inerts from decay (g VSS/g VSS)
% sedim: sedimentation coefficient (-)
% Y_DO_a_g: Utilization rate of DO by algae during growth (g O2/g VSS)
% Y_DO_a_d: Utilization rate of DO by algae during decay (g O2/g VSS)
% Y_IC_a_g: Utilization rate of IC by algae during growth (g C/g VSS)
% Y_IC_a_d: Utilization rate of IC by algae during decay (g C/g VSS)
% Y_N_a_g: Utilization rate of AIN by algae during growth (g N/g VSS)
% Y_N_a_d: Utilization rate of AIN by algae during decay (g N/g VSS)

    % Bacteria associated

% Y_20 : Yield bacterial growth at at 20*C (g VSS/g bCOD)   
% mu_20 : maximal specific growth at 20*C (g VSS/g VSS/s)
% Ks_20 : affinity on bCOD at 20*C (mg bCOD/L)
% kd_20 : decay rate (g VSS/g VSS/s)
% teta_mu : correction coefficient of max specific growth for temperature
% teta_Ks : correction coefficient of affinity on bCOD for temperature
% teta_kd : correction coefficient of decay rate for temperature
% K_DO_bacteria : affinity on dissolved oxygen(mg/L)
% NO3: NO3 concentration due to nitrifier activity (mg NO3/L)
% Y_DO_b_g : dissolved oxygen consuption associated to bacterial growth (g O2/gVSS)
% Y_C_b_g : inorganic carbon consumption associated to bacterial growth (g C/gVSS)
% Y_DO_b_d : dissolved oxygen consumed associated to bacterial decay (g O2/g VSS)
% Y_C_b_d : Inorganic carbon production associated to bacterial decay (g C-CO2/g VSS)
% Y_N_b_g : Utilization rate of AIN by heterotrophic bacteria during growth (g N/g VSS)
% Y_N_b_d : Utilization rate of AIN by heterotrophic bacteria during decay (g N/g VSS)
% Y_DO_n: DO consumption from nitrifiers activity (g DO/g VSS)
% Y_C_n: Inorganic carbon consumption from nitrifiers activity (g C/g VSS)
% Y_alk_n: inert charges consumption from nitrifiers activity (mol eq /g VSS)

    % Equilibrium constants

% Kla_O2: mass transfer coef of oxygen from the atm to the pond (s-1)
% Kla_CO2 : mass transfer coef of CO2 form atmosphere to water (s-1)
% KP1 : equilibrium constant H3PO4 - H2PO4 (source = wiki)
% KP2 : equilibrium constant H2PO4 - HPO4
% KP3: equilibrium constant HPO4 - PO4
% KN: equilibrium constant NH4 - NH3

    % Inlet characteristics

% COD_IN: influent chemical oxygen demand (mg/L)
% IC_IN : Total inorganic carbon in influent (mg C/L)
% IP_IN : Total inorganic phosphorous in influent (mg P/L)
% IN_IN: Total amoniacal nitrogen in influent (mg N/L)
% Alk: inert charges balance in the inlet (mol eq/m3)
% T_inlet : Temeprature of inlet (C)
% X_bacteria_IN : Heterotrophic bacteria concentration in the inlet (mg/L)
% C_IN : concentration of influent in total coliforms supposed constant(MPN/m3)

    % Miscellaneous
    
% sigma_a: transmitance of algal broth (m-1)
% alpha_dis: proportionality factor for disinfection rate against sunlight intensity (d-1.m2.J-1)                                                                 

    % Design parameters

% s_t: time difference between the input time and the output time (s)
% l_x : width of the pond (m)
% l_y : length of the pond (m)
% d_i : depth of the pond (m) 
% d_obj: design depth of the pond (m)
% HRT: Hydraulic Retention Time (d)
% eps_i: factor for inlet flowrate in case of non-continuous operations
% eps_o: factor for outlet flowrate in case of non-continuous operations
% option_vol_corr




%% OUTPUTS:

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


%% PRIMO CALCULATIONS

% S: pond surface (m2)
% V: pond volume (m3)
% qi: inlet flow rate calculated based on instant volume and actual depth -potentially different than design depth (m3/s)
% bCOD_IN : influent biological chemical oxygen demand (mg/L)

S = l_x*l_y;
V = S*d_i;
qi = qi_i;
bCOD_IN = 0.9*COD_IN;
 
%% Calculation of pH

pH = pH_pond_4( IC_i , IP_i, IN_i, Sigma_i/1000, Tp_i, KP2, KP3, KN );  

% Determination of temperature dependent equilibrium constants

KC1 = K_carbonate_1(Tp_i);
KC2 = K_carbonate_2(Tp_i);
Ke = K_ionic_product(Tp_i);

% Calculation of species concentration

nh = 10^(-pH);
H = nh*1000;                                                                              % conversion for mol-H/L to g-H/m3
OH = Ke/nh*1000;

H2CO3 = IC_i/(1 + KC1/nh + KC1*KC2/nh^2);
HCO3 = H2CO3*KC1/nh;
CO3 = H2CO3*KC1*KC2/(nh^2);

H2PO4 = IP_i/(1 + KP2/nh + KP2*KP3/nh);
HPO4 = H2PO4*KP2/nh;
PO4 = H2PO4*KP2*KP3/(nh^2);

NH4 = IN_i/(1 + KN/nh);
NH3 = NH4*KN/nh;

%% Pond Temperature determination

[Tp, me, T_soil,~,~,~,~,~,~,~ ] = pond_temperature_instant(Tp_i,Ta_i,Hs_dir_i,Hs_dif_i,RH_i,qr_i,Dw_a,Nu,T_soil_i,T_inlet,qi,s_t,V,S,l_y,rho_w,Cp_w,eps_w,sigma_stephan,fa,eps_a,Nu_a,R,Mw,Lw,lambda_a,alpha_a,alpha_s,Ts_ref,ks);
                  
%% Volume variation
% qe: water gain/loss due to condensation/evaporation
% diff: water gain/loss through inlet, evaporation, and rain (m3)
% DF: dilution factor due to rain and evaporation


qe = me/rho_w*S;
qi = eps_i*qi;

if d_i < 0.8*d_obj                                                                      %if the depth falls below a limit, a volume of wastewater is added for compensation
    V_add = (d_obj - d_i) * S;                                                    % leave the possibility here to take a volume added as an ouput and compute the water footprint of the algal pond.
else
    V_add = 0;
end
qi = qi + V_add/s_t;                                                                  % V_add is a volume added of wastewater so it is added to qi


if d_i > 2*d_obj                                                                         % correction if depth is getting to high only due to qi and qe
    if option_vol_corr == 1
        V_out = (d_i - d_obj)*S;
    else
        V_out = 0;                                                                          % No correction if option _vol_corr is  0
    end
else
    V_out = 0;
end

if option_vol_corr == 0                                                              % volume out forced (i.e. not flowing naturally)
    V_out = 0;
else
    V_out = V_out + eps_o*S*d_obj*1/HRT_i;                        % Volume withdrawn in case of manual operation (option_vol_corr != 0)
end

diff = (qi + qr_i - qe)*s_t - V_out; 
DF = (d_i + (qr_i - qe)*s_t/S)/d_i; 


% Adjustment of depth and calculation of qo

if option_vol_corr == 0
    if diff > 0
        if d_i < d_obj
            d = min(d_obj , d_i + diff/S);
            qo = max( 0 , (d_i + diff/S - d_obj)*S/s_t);
        else
            d = d_obj;
            qo = qi + qr_i - qe;
        end
    else
        d = d_i + diff/S;
        qo = 0;
    end
else
    d = d_i + diff/S;
    qo = V_out/s_t;
end
    
    
%% Bilogical variables        
    % Calculation of the evolution of concentration of the different
    % micro-organism (and bCOD) by mass-balance + growth kinetics
    
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
    % determination of the concentration evolution of DO, IC, IN, IP,
    % charge balance accounting for biological assimilation/utilization,
    % gas transfer, mass balance.

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
[kd , k_nat , k_pH , k_sun] = PhD_disinfection(s_t , Tp_i , pH , G_0 , sigma_a , d_i , alpha_dis); % instant decay rate (s-1)
Kd = kd/24/3600;


C_coli = C_coli_i/DF - Kd*C_coli_i*s_t + (C_IN - sedim*C_coli_i)*qi*s_t/V;

E_coli_kill_nat = k_nat*C_coli_i*s_t*V;
E_coli_kill_pH = k_pH*C_coli_i*s_t*V;
E_coli_kill_sun = k_sun*C_coli_i*s_t*V;



end




