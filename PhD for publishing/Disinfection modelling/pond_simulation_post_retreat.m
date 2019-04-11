function [ Time,Hs_dir,Hs_dif,Tp,s_t,me,qr,X_algae,C_coli,DO,bCOD,X_bacteria,H2CO3,HCO3,CO3,TC,H,OH,H3PO4,H2PO4,HPO4,PO4,TP,pH,Gamma,X_algae_debris,X_bacteria_debris,X_nit,Carbon_IN,Carbon_stock,Carbon_OUT,Carbon_acc,O2_IN,O2_stock,O2_OUT,Gamma_2 ] = pond_simulation_post_retreat( str,option_dis,option_pH,p_t,n_ref,l_x,l_y,d,HRT,rho_w,Cp_w,eps_w,sigma_stephan,eps_a,fa,Lw,Nu_a,R,Mw,lambda_a,alpha_a,Cp_s,rho_s,ks,r,Kla_O2,Kla_CO2,Alk,KP1,KP2,KP3,K_C,K_DO_algae,sigma,PE,HV,eps_eff,lambda_algae_light,lambda_algae_dark,Y_DO_a_g,Y_DO_a_d,Y_IC_a_g,Y_IC_a_d,Y_alk_a_g,mu_20,Ks_20,Y_20,kd_20,teta_mu,teta_kd,teta_Ks,fs_i,K_DO_bacteria,Y_DO_growth,Y_C_growth,Y_DO_decay,Y_C_decay,Y_alk_growth,T_IN,COD_IN,C_IN,TC_IN,TP_IN,X_bacteria_IN,Tp_i,X_i,C_coli_i,DO_i,pH_i,TC_i,TP_i,bCOD_i,X_bacteria_i,Alk_i,k_s_m,k_pH_m,k_20_m,teta_m,k_DO_m, sedim ,NO3, Y_O2_nit, Y_C_nit,Y_alk_nit,X_nit_i, sigma_a,alpha_dis )

% This function aims at performing a full simulation of an HRAP environmental parameters
% and disinfection performance from design parameters, influent characteristics, and meteorological
% data.

%% INPUTS

    % Meteorological conditions:

% xls file taken from a string (argument str)

    % Design parameters

% l_x : width of the pond (m)
% l_y : length of the pond (m)
% d : depth of the pond (m) 
% HRT : Hydraulic Retention Time (d)

    % Calculation options:

% option_dis : choice for disinfection modelling: 1 = Marais (1974); 2 = Mayo
%   (1995) ; 3 = Craggs (2003) ; 4 = manual adjustment; 5 = PhD
%   disinfection
% option_pH: choice for the calculation of pH: 1 = the pH is not included in the simulation 
%   (technically it is assuming it is constant = 7) and inorganic carbon is tweaked to be in excess 
%   (play on TCin); 2 = the pH is calculated accounting only for carbonate
%   species equilibria; 3 = the pH is calcualted accounting for carbonate and phosphorous species equilibria
% p_t : factor division for the step of time in pond temperature calcualtions (being one hour if p_t = 1)
% n_ref : division factor for the refinement of calculations in case of negative inorganic carbon 

    % Physico-chemical constants

% rho_w : density of pond water (kg/m3)
% Cp_w : speciffic heat capacity of pond water (J/kg/K)
% V : pond volume (m3)
% eps_w : water emissivity
% sigma_stephan : Stephan-Boltzmann constant (W/m�/K) 
% S : pond surface (m�)
% fa : algal absorption fraction (%)
% eps_a : air emissivity
% Lw : water latent heat (J/kg)
% Nu_a : air kinematic viscosity (m�/s)
% R : ideal gas constant (Pa.m3/mol/K)
% Mw : molecular weight of water (kg/mol)
% lambda_a : air thermal conductivity (W/m/k)
% alpha_a : air thermal diffusivity (m�/s)
% Cp_s : soil specific heat capacity ((J/kg/K)
% rho_s : soil density (kg/m3)
% ks : soil thermal conductivity (W/m/K)
% r : reflexivity of the pond (%)
% Kla_O2: mass transfer coef of oxygen from the atm to the pond (s-1)
% Kla_CO2 : mass transfer coef of CO2 form atmosphere to water (s-1)
% H2CO3_d_sat : saturation concentration of dissolved carbon dioxide (mg C-CO2/L)
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
% lambda_algae: rate of respiration of algae (kg TSS decayed/kg TSS produced/s). 1 value in the dark, one value in the light; 
% Y_DO_a_g: Yield of DO production associated with algal growth (g DO/g TSS produced)
% Y_DO_a_d: Yield of DO consumption associated with algal decay (g DO/g TSS decayed)
% Y_IC_a_d: Yield of inorganic carbon consumption associated with algal growth (g C/g TSS decayed)
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
% K_DO : affinity on dissolved oxygen(mg/L)
% Y_DO_growth : dissolved oxygen consuption associated to bacterial growth
% (g O2/gVSS)
% Y_HCO3_growth : bicarbonated consumption associated to bacterial growth
% (g C-HCO3-/gVSS)
% Y_CO2_growth : Carbon dioxide consumption associated to bacterial growth
% (g C-CO2/g VSS)
% Y_DO_decay : dissolved oxygen consumed associated to bacterial decay (g
% O2/g VSS)
% Y_CO2_decay : Carbon dioxide production associated to bacterial decay (g
% C-CO2/g VSS)
% Y_alk_growth : alkalinity contribution through bacterial growth (mol eq/g
% VSS)

% NO3 : NO3 concentration due to nitrifier activity (mg NO3/L)

    % Inlet characteristics

% T_IN : Temeprature of inlet (C)
% COD_IN: influent chemical oxygen demand (mg/L)
% C_IN : concentration of influent in total coliforms supposed constant(MPN/m3)
% TC_IN : Total inorganic carbon in influent (mg C/L)
% TP_IN : Total inorganic phosphorous in influent (mg P/L)
% X_bacteria_IN : Heterotrophic bacteria concentration in the inlet (mg/L)


    % Initialisation

% Tp_i : Initial temperature of the pond (C)
% X_i : initial concentration in algae (kg/m3)
% C_coli_i : initial concentration in E. coli (MPN/m3)
% DO_i : initial dioxygen concentration (mg/L)
% pH_i : initial pH
% TC_i : inital total inorganic carbon concentration (mg/L
% TP_i : initial total phopshorous (mg/L)
% bCOD_i : initial bCOD concentration (mg/L)
% X_bacteria_i : initial heterotrophic bacteria concentration (g/m3)
% Alk_o : initial pond alkalinity (mol/L)

    % Manual disinfection rate
    
% k_s_m: solar dependent manually adjusted coefficient (m2.W-1.d-1)
% k_pH: pH dependent manually adjusted coefficient (d-1)
% k_20_m: temperature dependent manually adjusted coefficient (d-1)
% teta_m: temperature dependent Arrheniusish manually adjusted coefficient (-)
% k_DO_m: Disolved oxygen concentration dependent manually adjusted coefficient (L.mgDO-1.d-1)

%% OUTPUTS

% Time : Table of time of the year (s)
% Hs_dir : Table of direct solar radiation or global if assuming it is only direct according to Time (W/m2)
% Hs_dir : Table of diffuse solar radiation according to Time (W/m2)
% Tp : Table of temperature of the pond according to Time (K)
% s_t : step of time (s)
% me : mass of evaporation (kg/m2/s)
% qr : rain water flow (m3/m2/s)
% X : algal concentration (kg/m3)
% C_coli: concentration in total colifoms (MPN/m3)
% DO: oxygen concentration in pong (mg/L)
% bCOD: biological chemical oxygen demand in the pond (mg/L) taken as substrate indicator for
% heterotrophic bacteria
% X_bacteria : heterotrophic bacteria population (mg/L)
% KC1 : dissiciation constant for acido-basic reaction between dissolved CO2
% and bicarbonate (mol/L)
% KC2 : dissociation constant for acido-basic reaction between bicarbonate
% and carbonate (mol/L)
% Ke : dissociation constant for water ionic product (mol^2/L^2)
% H2CO3 : concentration of aqueous carbon dioxide in the pond (mg C-CO2/L)
% HCO3 : concentration of bicarbonate in the pond (mg C-HCO3-/L)
% CO3 : concentration of carbonate in the pond (mg C-CO32-/L)
% H : concentration of H+ proton in the pond (mol/m3)
% OH : concentration of OH- in the pond (mol/m3)
% H3PO4 : concentration of H3PO4 in the pond (mg P-H3PO4/L)
% H2PO4 : concentration of H2PO4- in the pond (mg P-H2PO4-/L)
% HPO4 : concentration of H(PO4)2- in the pond (mg P-H(PO4)2-/L)
% PO4 : concentration of (PO4)3- in the pond (mg P-(PO4)3-/L)
% TP : total phosphorous concentration (mg P/L)
% pH : pH in the pond (dimensionless)
% Gamma : Alkaline / ionic balance of the broth (meq/L)
% X_debris: biological debris inert concentration (mg/L)
% X_nit: nitrifying bacteria (mg/L)
% Carbon_IN: mass of carbon entering the pond (g)
% Carbon_stock: mass of carbon inside the pond (g)
% Carbon_OUT: mass of carbon leaving the pond (g)

%% Function


% PRIMO CALCULATIONS

% S: pond surface (m2)
% V: pond volume (m3)
% qi: inlet flow rate (m3/s)
% qo: outlet flow rate (m3/s)
% bCOD_IN : influent biological chemical oxygen demand (mg/L)

S = l_x*l_y;
V = S*d;
qi = V/HRT/24/3600;
qo = qi;
bCOD_IN = 0.9*COD_IN;

% Import data

xlsread(char(str));



% Temperature calculation

[Time, Hs_dir, Hs_dif, Tp, s_t, me, qr] = pond_temperature(Tp_i,T_IN,rho_w,Cp_w,V,eps_w,sigma_stephan,S,fa,eps_a,Lw,l_y,Nu_a,R,Mw,lambda_a,alpha_a,Cp_s,rho_s,ks,qi,strcat(str,'_model'),p_t);

% Creation of the rest of the outputs

n = size(Time);
n = n(2);

X_algae = zeros(n,1);
C_coli = zeros(1,n);
DO = zeros(1,n); 
pH = zeros(1,n); 
H2CO3 = zeros(1,n); 
HCO3 = zeros(1,n); 
CO3 = zeros(1,n); 
TC = zeros(1,n);
OH = zeros(1,n);
H = zeros(1,n);
bCOD = zeros(1,n);
X_bacteria = zeros(1,n);
H3PO4 = zeros(1,n);
H2PO4 = zeros(1,n);
HPO4 = zeros(1,n);
PO4 = zeros(1,n);
TP = zeros(1,n);
Gamma = zeros(1,n);
X_algae_debris = zeros(1,n);
X_bacteria_debris = zeros(1,n);
X_nit =  zeros(1,n);
Carbon_IN = zeros(1,n);
Carbon_stock = zeros(1,n);
Carbon_OUT = zeros(1,n);
Carbon_acc = zeros(1,n);
O2_IN = zeros(1,n);
O2_stock = zeros(1,n);
O2_OUT = zeros(1,n);
Gamma_2 = zeros(1,n);

% Initialisation:

X_algae(1) = X_i;
C_coli(1,1) = C_coli_i;
KC1 = K_carbonate_1(Tp(1));
KC2 = K_carbonate_2(Tp(1));
Ke = K_ionic_product(Tp(1));
DO(1) = DO_i;
pH(1) = pH_i;
H(1) = 10^(-pH(1))*1000;
OH(1) = 10^6*Ke/H(1);
TC(1) = TC_i;
H2CO3(1) = TC(1)/( 1 + KC1/(H(1)/1000) + KC1*KC2/(H(1)/1000)^2 );
HCO3(1) = 1000*KC1*H2CO3(1)/(H(1));
CO3(1) = 1000*KC2*HCO3(1)/(H(1));
bCOD(1) = bCOD_i;
X_bacteria(1) = X_bacteria_i;
TP(1) = TP_i;
H3PO4(1) = TP(1)/(1 + KP1/(H(1)/1000) + KP1*KP2/(H(1)/1000)^2 + KP1*KP2*KP3/(H(1)/1000)^3);
H2PO4(1) = 1000*H3PO4(1)*KP1/H(1);
HPO4(1) = 1000*H2PO4(1)*KP2/H(1);
PO4(1) = 1000*HPO4(1)*KP3/H(1);
Gamma(1) = Alk_i;
Gamma_2(1) = Alk_i;
X_nit(1) = NO3/sedim*0.0196/0.98*113/62;

% Disinfection rate initialisation:

if option_dis <= 4
    [ k_s, k_pH,k_20,teta,k_DO ] = HRAP_disinfection_rate( option_dis,d,k_s_m,k_pH_m,k_20_m,teta_m,k_DO_m );
    Kd = 0;
else
    Kd = PhD_disinfection(Time(2)-Time(1),Tp(1),pH(1),Hs_dir(1) + Hs_dif(1),sigma_a,d,alpha_dis);
    k_s = 0;
    k_pH = 0;
    k_20 = 0;
    teta = 1;
    k_DO = 0;
end
    
i = 1;
alpha = [1,0]

a_1 = 0;
a_2 = 0;
a_3 = 0;
a_4 = 0;
a_5 = 0;
a_6 = 0;
a_7 = 0;

while i <= n-1
    
    refine_loop = 1;
    
    if TC(i) < 0 || DO(i) < 0
        
        refine_loop = refine_loop + 1;
        if refine_loop > 5
            break
        end
        
        i = i-1;
        
        Time = refinement_linear_interpol(Time,n_ref,i);
        Hs_dir = refinement_linear_interpol(Hs_dir,n_ref,i);
        Hs_dif = refinement_linear_interpol(Hs_dif,n_ref,i);
        me = refinement_linear_interpol(me,n_ref,i);
        Tp = refinement_linear_interpol(Tp,n_ref,i);
%         C_IN_table = refinement_linear_interpol(C_IN_table,n_ref,i);
%         qr = refinement_linear_interpol(qr,n_ref,i);
        
        Add = zeros(1,n_ref -1);
        
        DO = [DO,Add];
        pH = [pH,Add];
        H2CO3 = [H2CO3,Add];
        HCO3 = [HCO3,Add];
        CO3 = [CO3,Add];
        TC = [TC,Add];
        OH = [OH,Add];
        H = [H,Add];
        H3PO4 = [H3PO4,Add];
        H2PO4 = [H2PO4,Add];
        HPO4 = [HPO4,Add];
        PO4 = [PO4,Add];
        TP = [TP,Add];
        bCOD = [bCOD,Add];
        X_bacteria = [X_bacteria,Add];
        X_algae = [X_algae; Add'];
        C_coli = [C_coli,Add];
        Gamma = [Gamma, Add];
        X_algae_debris = [X_algae_debris, Add];
        X_bacteria_debris = [X_bacteria_debris, Add];
        X_nit = [X_nit, Add];
        Carbon_IN = [Carbon_IN, Add];
        Carbon_stock = [Carbon_stock, Add];
        Carbon_OUT = [Carbon_OUT, Add];
        Carbon_acc = [Carbon_acc,Add];
        O2_IN = [O2_IN, Add];
        O2_stock = [O2_stock, Add];
        O2_OUT = [O2_OUT, Add];
        Gamma_2 = [Gamma_2, Add];
        n = n + n_ref - 1;
    end
    
    s_t = Time(i+1) - Time(i);
    
    
    % Algae concentration evolution
    [ P_algae_gross, P_algae_resp, P_algae_net ] = algae_simple_productivity( Hs_dir(i), Hs_dif(i), sigma, X_algae(i), d, S, TC(i), K_C, PE, HV, eps_eff,DO(i),K_DO_algae,lambda_algae_light,lambda_algae_dark);

    % Change for Benoit of respiration (01/09/2017)
    kd_a = kd_20*teta_kd^(Tp(i) - 273.15 - 20);
    P_algae_resp = (1 - fs_i)*(DO(i)/(K_DO_algae + DO(i)))*kd_a*X_algae(i)*d*S;
    Q_algae_debris = fs_i*(DO(i)/(K_DO_algae + DO(i)))*kd_a*X_algae(i)*d*S;
    
    P_algae_net = P_algae_gross - P_algae_resp - Q_algae_debris;
    
    X_algae(i+1) = X_algae(i) + P_algae_net*s_t/V - sedim*qo*X_algae(i)*s_t/V;
    
    % Bacteria concentration evolution
    [ Q_bacteria_gross, Q_bacteria_decay ] = heterotrophic_bacteria_growth( Tp(i),mu_20,Ks_20,kd_20,teta_mu,teta_Ks,teta_kd,DO(i),K_DO_bacteria,bCOD(i),X_bacteria(i),d,S );
    Q_bacteria_debris = fs_i*Q_bacteria_decay;
    Q_bacteria_decay = (1-fs_i)*Q_bacteria_decay;    
    X_bacteria(i+1) = X_bacteria(i) + qi*s_t*X_bacteria_IN/V - sedim*qo*s_t*X_bacteria(i)/V + Q_bacteria_gross*s_t/V - Q_bacteria_decay*s_t/V - Q_bacteria_debris*s_t/V;

    % bCOD concentration evolution
    bCOD(i+1) = bCOD(i) - (1/Y_20)*Q_bacteria_gross*s_t/V + qi*bCOD_IN/V*s_t - sedim*qo*bCOD(i)*s_t/V ;    %     bCOD(i+1) = bCOD(i) - (1/Y_20)*Q_bacteria_gross*s_t/V + qi*bCOD_IN/V*s_t - sedim*qo*bCOD(i)*s_t/V + (1-fs_i)/Y_20*Q_bacteria_decay*s_t/V;

    % Debris concentration evolution (mg VSS/L)

    X_algae_debris(i+1) = X_algae_debris(i) + Q_algae_debris*s_t/V - sedim*qo*X_algae_debris(i)*s_t/V; %gTSS/L
    X_bacteria_debris(i+1) = X_bacteria_debris(i) + Q_bacteria_debris*s_t/V - sedim*qo*X_bacteria_debris(i)*s_t/V;
    X_debris(i+1) = X_algae_debris(i+1) + X_bacteria_debris(i+1);
    
    % Nitrifiers
    P_NO3 = NO3*qo; % NO3 productivity (g NO3/s)
    Q_nit = P_NO3*0.0196/0.98*113/62; % g VSS nit/s
    X_nit(i+1) = X_nit(i) + Q_nit*s_t/V - sedim*qo*X_nit(i)*s_t/V;%neglecting nitrifiers debris and decay
    
    % Calculation of DO concentration (mg/L)
    O2_algae_net = 1000*(Y_DO_a_g*P_algae_gross - Y_DO_a_d*P_algae_resp); % g-O2/s    
    DO(i+1) = DO(i) + s_t*(saturation_DO_no_salinity(Tp(i)) - DO(i))*Kla_O2 + O2_algae_net*s_t/V - qo*s_t*DO(i)/V + 0*qi*s_t*saturation_DO_no_salinity(Tp(i))/V - Y_DO_growth*Q_bacteria_gross*s_t/V - Y_O2_nit*Q_nit*s_t/V - Y_DO_decay*Q_bacteria_decay*s_t/V; % Q_bacteria_gross*Y_DO_growth*s_t/V - Q_bacteria_decay*Y_DO_decay*s_t/V ;
    
    %     DO(i+1) = DO(i) + s_t*(saturation_DO_no_salinity(Tp(i)) - DO(i))*Kla_O2 + O2_algae_net*s_t/V - qo*s_t*DO(i)/V + 0*qi*s_t*saturation_DO_no_salinity(Tp(i))/V + 1.42*(Q_bacteria_gross - Q_bacteria_decay)*s_t/V - (bCOD_IN*qi - bCOD(i)*qo)/V*s_t - Q_O2_nit*s_t/V; % Q_bacteria_gross*Y_DO_growth*s_t/V - Q_bacteria_decay*Y_DO_decay*s_t/V ;
    
    % Calculations of inorganic TC and TP to prepare the calcualtion of pH:
    H2CO3_d_sat = CO2_solubility(Tp(i));
    TC_supply = TC_IN*qi*s_t/V + Kla_CO2*(H2CO3_d_sat - H2CO3(i))*s_t - Y_IC_a_g*P_algae_gross*1000*s_t/V + Y_IC_a_d*P_algae_resp*1000*s_t/V - qo*H2CO3(i)*s_t/V  - qo*s_t*HCO3(i)/V - qo*s_t*CO3(i)/V + Y_C_growth*Q_bacteria_gross*s_t/V - Y_C_nit*Q_nit*s_t/V + Y_C_decay*Q_bacteria_decay*s_t/V; %term from carbon realeased by biomass creation
    %- Q_bacteria_gross*Y_HCO3_growth*s_t/V ; %+ Q_bacteria_gross*Y_CO2_growth*s_t/V + Q_bacteria_decay*Y_CO2_decay*s_t/V 
    %     TC_supply =  TC_IN*qi*s_t/V + Kla_CO2*(H2CO3_d_sat - H2CO3(i))*s_t - Y_IC_a_g*P_algae_gross*1000*s_t/V + Y_IC_a_d*P_algae_resp*1000*s_t/V - qo*H2CO3(i)*s_t/V  - qo*s_t*HCO3(i)/V - 12/5/8*(1.42*(Q_bacteria_gross - Q_bacteria_decay)*s_t/V - (bCOD_IN*qi - bCOD(i)*qo)/V*s_t) - 1/113*20*(Q_bacteria_gross - Q_bacteria_decay)*2.5*12/50*s_t/V - qo*s_t*CO3(i)/V; %- Q_bacteria_gross*Y_HCO3_growth*s_t/V ; %+ Q_bacteria_gross*Y_CO2_growth*s_t/V + Q_bacteria_decay*Y_CO2_decay*s_t/V 
    TC(i+1) = TC(i) + TC_supply;
    
    P_supply = -0.01*(P_algae_gross - P_algae_resp)*s_t*1000/V - 0.01*(Q_bacteria_gross - Q_bacteria_decay + Q_nit)*s_t/V + TP_IN*qi*s_t/V - TP(i)*qo*s_t/V; 
    TP(i+1) = TP(i) + P_supply;
    
    % Calculation of alkaline balance evolution
    
    Gamma(i+1) = Gamma(i) + qi*s_t*Alk/V - qo*s_t*Gamma(i)/V + Y_alk_a_g*P_algae_gross*s_t*1000/V + Y_alk_growth*Q_bacteria_gross*s_t/V + Y_alk_nit*Q_nit*s_t/V;
    Gamma_2(i) = (HCO3(i) + 2*CO3(i))/12*1000/1000 + OH(i)/1000*1000 - H(i)/1000*1000;
    
    % pH calculation
    pH(i+1) = pH_pond( option_pH, TC(i+1) , TP(i+1), Gamma(i+1)/1000, Tp(i+1), KP1, KP2, KP3 );

    % Concentration of species in equilibrium calculations
    KC1 = K_carbonate_1(Tp(i+1));
    KC2 = K_carbonate_2(Tp(i+1));
    Ke = K_ionic_product(Tp(i+1));
    
    nh = 10^(-pH(i+1));
    H(i+1) = nh*1000; % conversion for mol-H/L to g-H/m3
    OH(i+1) = Ke/nh*1000;

    H2CO3(i+1) = TC(i+1)/(1 + KC1/nh + KC1*KC2/nh^2);
    HCO3(i+1) = H2CO3(i+1)*KC1/nh;
    CO3(i+1) = H2CO3(i+1)*KC1*KC2/(nh^2);

    H3PO4(i+1) = TP(i+1)/(1 + KP1/nh + KP1*KP2/nh^2 + KP1*KP2*KP3/nh^3);
    H2PO4(i+1) = H3PO4(i+1)*KP1/nh;
    HPO4(i+1) = H3PO4(i+1)*KP1*KP2/(nh^2);
    PO4(i+1) = H3PO4(i+1)*KP1*KP2*KP3/(nh^3);
    
    % E. coli disinfection calculations
    
    G_0 = (Hs_dir(i) + Hs_dif(i));
    if option_dis <= 4
        Kd = (k_s*G_0/d + k_pH*pH(i) + k_20*teta^(Tp(i) - 273.15 - 20) + k_DO * DO(i))/24/3600; % instant decay rate (s-1)
    else
        Kd = PhD_disinfection(Time(i+1)-Time(i) , Tp(i) , pH(i) , G_0 , sigma_a , d , alpha_dis)/24/3600; % instant decay rate (s-1)
        Kd*24*3600
    end
        
    C_coli(i+1) = C_coli(i) - Kd*C_coli(i)*s_t + (C_IN - sedim*C_coli(i))*qi*s_t/V;
    
%     %% Balance calculations
%     % Carbon (g)
% 
% 
%     Carbon_IN(i) = qi*TC_IN*s_t + max(0,Kla_CO2*(H2CO3_d_sat - H2CO3(i))*s_t*V) + 0.300*bCOD_IN*qi*s_t + 0.5310*X_bacteria_IN*qi*s_t;
%     Carbon_stock(i) = (H2CO3(i) + HCO3(i) + CO3(i))*V + 0.5096*(X_algae(i) + X_algae_debris(i))*V*1000 + 0.5310*(X_bacteria(i) + X_bacteria_debris(i) + X_nit(i))*V + 0.300*bCOD(i)*V ;
%     Carbon_OUT(i) = max(0,-Kla_CO2*(H2CO3_d_sat - H2CO3(i))*s_t*V) + (H2CO3(i) + HCO3(i) + CO3(i))*qo*s_t + sedim*(0.5096*(X_algae(i) + X_algae_debris(i))*1000 + 0.5310*(X_bacteria(i) + X_bacteria_debris(i) + X_nit(i)) + 0.300*bCOD(i))*qo*s_t ;
%     
%     if i>= 2
%         a_7 = a_7 + Carbon_stock(i) - Carbon_stock(i-1) - Carbon_IN(i-1) + Carbon_OUT(i-1)
%     end
%     Carbon_acc(i) = Carbon_IN(i) - Carbon_OUT(i);
% 
%        
%     a_1 = a_1 + X_algae(i+1)*V - X_algae(i)*V - (P_algae_net*s_t - sedim*qo*s_t*X_algae(i))
%     a_2 = a_2 + X_bacteria(i+1)*V - X_bacteria(i)*V - (Q_bacteria_gross*s_t - Q_bacteria_decay*s_t - Q_bacteria_debris*s_t + X_bacteria_IN*qi*s_t - sedim*X_bacteria(i)*qo*s_t) 
%     a_3 = a_3 + X_debris(i+1)*V - X_debris(i)*V - ((Q_algae_debris + Q_bacteria_debris)*s_t - sedim*X_debris(i)*s_t*qo)
%     a_4 = a_4 + X_nit(i+1)*V - X_nit(i)*V - (Q_nit*s_t - sedim*X_nit(i)*qo*s_t)
%     a_5 = a_5 + bCOD(i+1)*V - bCOD(i)*V - (  qi*bCOD_IN*s_t  - (1/Y_20)*Q_bacteria_gross*s_t - sedim*qo*bCOD(i)*s_t)
    
    % Oxygen balance:
    
%     O2_stock(i) = 0.3724*(X_algae(i)  + X_algae_debris(i))*V*1000 + 0.2832*(X_bacteria(i) + X_nit(i) + X_bacteria_debris(i))*V + 0.1200*bCOD(i)*V + DO(i)*V;
%     O2_IN(i) = 0*qi*s_t*saturation_DO_no_salinity(Tp(i)) + max(0,s_t*(saturation_DO_no_salinity(Tp(i)) - DO(i))*Kla_O2*V) + 0.1200*bCOD_IN*qi*s_t + 0.2832*X_bacteria_IN*qi*s_t;
%     O2_OUT(i) = max(0,-s_t*(saturation_DO_no_salinity(Tp(i)) - DO(i))*Kla_O2*V) + sedim*(0.3724*(X_algae(i)  + X_algae_debris(i))*1000 + 0.2832*(X_bacteria(i) + X_nit(i) + X_bacteria_debris(i)) + 0.1200*bCOD(i))*qo*s_t + DO(i)*qo*s_t ;
%     
%     if i>=2
%         a_6 = a_6 + O2_stock(i) - O2_stock(i-1) + O2_IN(i-1) - O2_OUT(i-1)
%     end
    
    %% end of lo.op
    
    i = i+1;
%     alpha(2) = (Time(end) - Time(i))/24/3600
    
%     alpha(2) = (n-i)/n*100
    
    if alpha(2) == alpha(1)
        gamma = 1;
    end
    
    alpha(1) = alpha(2);
    

end

X_algae = X_algae';

end




