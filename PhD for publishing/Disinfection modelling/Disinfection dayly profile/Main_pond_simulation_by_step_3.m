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

p_t = 100;
 
f_refine = 3;
t_cut = 1; % step of time under which the computing is cut (s)

% When comparing with real data, base case is p_t = 100, f_refine = 3,
% t_cut = 0.01;

%% INPUTS

% Examples of input strings: Test_10_Feb_repeat - Palmy_year_1 -
% Palmy_year_2 - Palmy_year_3 - Test_16_Mar_repeat - Test_30_Sep_repeat -
% Test_28_Oct_repeat - Palmy_validation_year_2

str = 'Palmy_validation_year_2';
Time_start = datenum(2016,07,22,10,0,0);

% IMPORTANT MODELLING PARAMETERS (i.e. calculation implication, design
% parameters of HRAP, unceratain variables to set up)

l_x = 1;
l_y = 3.42;
d_obj = 0.25;

sedim = 0.8; % experimentally determined
PE = 0.02; % calibration
NO3 = 58;% mg/L, based on medium value

% Initialisation

HRT_1 = 7.9; % Median value over the year
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

%% Variables initialisation

% Meteo data
str = strcat(str,'_model');
raw_data = xlsread(char(str));

% n : step of time for calculation
% p : step of depth for soil temperature calculation
% s_t : step of time (s)
% p_t: initial interpolation of metorological data

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


for z = 1:n-1
    for k = 1:p_t
        A(p_t*(z-1) + k) = qr(z);
    end
end
qr = A;

n = (n-1)*p_t + 1;
s_t = s_t/p_t;
s_t_reference = s_t;


%% Inputs variable by time step

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

% Initialisation

Time(1) = Time_start;

% Design

d(1) = d_obj;
HRT(1) = HRT_1;

Kla_O2(1) = Kla_O2_1;
Kla_CO2(1) = Kla_O2(1)*(1.67/2.01)^(1/2); % source for diffusivity in water is engineeringtoolbox, values at 20°C
Kla_NH3(1) = Kla_O2(1)*(1.50/2.01)^(1/2);

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

Time(2) = Time(1) + s_t/24/3600;

refinement_index = 0;

i = 1;
while i < n

s_t = (Time(i+1) - Time(i))*24*3600;
if s_t < 0
    break
else
    if s_t < t_cut/f_refine
        break
    end
end
HRT(i) = HRT(1);
Kla_O2(i) = Kla_O2(1);
Kla_CO2(i) = Kla_CO2(1);
Kla_NH3(i) = Kla_NH3(1);

[ X_algae_c, bCOD_c , X_bacteria_c, C_coli_c  , X_debris_a_c , X_debris_b_c, X_nit_c , DO_c, pH_c, H2CO3_c , HCO3_c , CO3_c  , IC_c, H_c , OH_c , NH3_c, NH4_c, IN_c, H2PO4_c , HPO4_c , PO4_c , IP_c , Sigma_c  , kd_c , Tp_c , me_c , T_soil_c ,  d_c , qi_c , qo_c ] = ...
            pond_simulation_step_3( Hs_dir(i) , Hs_dif(i) , X_algae(i) , IC(i) ,DO(i) , bCOD(i) , X_bacteria(i) , IP(i) ,  Sigma(i) , Ta(i) , RH(i) , qr(i) , T_soil(i,:)', Tp(i), X_debris_a(i) , X_debris_b(i), X_nit(i), C_coli(i), IN(i) ...
            , Dw_a(i) , Nu(i) , rho_w , Cp_w , eps_w , sigma_stephan , fa , eps_a , Lw , Nu_a , R , Mw , lambda_a , alpha_a , alpha_s , Ts_ref , ks ...
            , K_C , PE , HV , eps_eff , K_DO_algae ,lambda_algae_light , lambda_algae_dark , sigma_a , fs_i , sedim ,  Y_DO_a_g , Y_DO_a_d , Y_C_a_g , Y_C_a_d , Y_N_a_g, Y_N_a_d  ...
            , Y_20 , mu_20 , Ks_20 , kd_20 , teta_mu , teta_Ks , teta_kd , K_DO_bacteria , NO3 , Y_DO_b_g ,  Y_DO_b_d , Y_C_b_g , Y_C_b_d , Y_N_b_g , Y_N_b_d , Y_DO_n, Y_C_n, Y_N_n, Y_alk_n ...
            , Kla_O2(i) , Kla_CO2(i) , Kla_NH3(i),  KP2, KP3 , KN ...
            , COD_IN , IC_IN , IP_IN , IN_IN, Sigma_IN , T_inlet, X_bacteria_IN , C_IN ...
            , sigma_a, alpha_dis ...
            , s_t , l_x , l_y , d(i), d_obj , HRT(i) );

    if i == 1
        if IC_c < 0 || DO_c < 0
            warning('/!\ Change initial conditions')
            break
        end
    end
    
    if IC_c < 0 || DO_c < 0

                      
        if (Time(i) - Time(i-1))*24*3600 > t_cut
            
            i = i-1;
            refinement_index = refinement_index + 1;
            
            s_t = (Time(i+1) - Time(i))*24*3600;
            Zeros = zeros( f_refine - 1 , 1);

            Aux = zeros( f_refine - 1 , 1);
            for k = 1:f_refine - 1
                Aux(k) = Time(i) + k*s_t/f_refine/24/3600;
            end

            Time = [ Time(1:i) ; Aux ; Time(i+1:end) ];



            Hs_dir = refinement_linear_interpol(Hs_dir,f_refine,i);
            Hs_dif = refinement_linear_interpol(Hs_dif,f_refine,i);
            Ta = refinement_linear_interpol(Ta,f_refine,i);
            RH = refinement_linear_interpol(RH,f_refine,i);
            Nu = refinement_linear_interpol(Nu,f_refine,i);
            Dw_a = -2.775*10^(-6)+4.479*10^(-8)*Ta + 1.656*10^(-10)*Ta.^2;
            qr = [qr(1:i) ; zeros(f_refine - 1,1) ; qr(i+1:end)];
            
            
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

            n = n + f_refine - 1;

            if Time(i+2) == 0
                Time(i+2) = Time(i+1) + s_t_reference/24/3600;
            end
        else
            if IC_c < 0
                IC_c = 0.001;
            end
            if DO_c < 0
                DO_c = 0.001;
            end
            
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

            if Time(i+2) == 0
                Time(i+2) = Time(i+1) + s_t_reference/24/3600;
            end

            i = i+1;
        end
    
    else
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

        if Time(i+2) == 0
            Time(i+2) = Time(i+1) + s_t_reference/24/3600;
        end

        i = i+1;
        if rem(i,10000) == 0
            i/n
        end
    end            
end

Time = Time(1:n);

%% Creation of estimates

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

t_Tp_20
t_pH_10
mean(max_Tp) - 273.15
mean(max_pH)

load gong.mat;
sound(y)

%%

%%


%% Comparing with real data
 
A1 = xlsread('\\tsclient\C\Users\pchambon\Dropbox\Pilot Scale\pH_DO_monitoring\HRAP A\HRAP_A_pH_DO_T_general_validation.xlsm','A2:F26639');

p = size(A1,1);
pH_meas = A1(1:p -1,2);
DO_meas = A1(1:p-1,4);
Tp_meas = A1(1:p-1,6);
Time_meas = A1(1:p-1,1);
Time_meas = Time_meas + 6.9396*10^5;

%% Figures

figure(100), clf, hold on
plot(Time,Tp-273.15,'k','LineWidth',1), ylabel({'Pond' 'Temperature'  '(^oC)'},'Fontsize',18); plot(Time_meas,Tp_meas,'ok');
ax = gca; xlim(gca,[Time(1) Time(end)]); datetick('x','dd-mmm-yyyy');
ax.FontSize = 14; 
ax.XTickLabelRotation = 45;
% datetick('x','dd-mmm-yyyy'); xlim(gca,[Time(850000) Time(end)]);

figure(101), clf, hold on
plot(Time,pH,'k','LineWidth',1), ylabel('pH','Fontsize',18);  plot(Time_meas,pH_meas,'ok');
ax = gca; datetick('x','dd-mmm-yyyy');
ax.FontSize = 14; ylim([6.0 12])
ax.XTickLabelRotation = 45;
% datetick('x','dd-mmm-yyyy'); xlim(gca,[Time(850000) Time(end)]);


figure(102), clf, hold on
plot(Time,DO,'k','LineWidth',1), ylabel({'Dissolved Oxygen'  '(mg.L^-^1)'},'Fontsize',18); plot(Time_meas,DO_meas,'ok');
ax = gca; xlim(gca,[Time(1) Time(end)]); datetick('x','dd-mmm-yyyy');
% ylim(gca,[0 18])
ax.FontSize = 14;
ax.XTickLabelRotation = 45;
% datetick('x','dd-mmm-yyyy'); xlim(gca,[Time(850000) Time(end)]);

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


% Figures over reduced periods
% n_1 = 527500;
% n_2 = n_1 + 30000;

% n_1 = 215000;
% n_2 = n_1 + 30000;

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
% datetick('x','dd-mmm-yyyy'); xlim(gca,[Time(850000) Time(end)]);

figure(122), clf, hold on
plot(Time,DO,'k','LineWidth',1), ylabel({'Dissolved Oxygen'  '(mg.L^-^1)'},'Fontsize',18); plot(Time_meas,DO_meas,'ok');
xlabel('July 2017','Fontsize',18)
ax = gca; datetick('x','dd-mmm-yyyy');
ax.FontSize = 14; 
ax.XTickLabelRotation = 45;
xlim(gca,[Time(n_1) Time(n_2)]);

% return

%% Y = X plot



Time_mod = [Time(1)];

Tp_mod_2 = [Tp(1)];
DO_mod_2 = [DO(1)];
pH_mod_2 = [pH(1)];

Time_meas_2 = [Time_meas(1)]
Tp_meas_2 = [Tp_meas(1)];
pH_meas_2 = [pH_meas(1)];
DO_meas_2 = [DO_meas(1)];

i = 1;
j = 2;
m = size(Time_meas,1);

while i <= n && j <= m
    if Time_meas(j) <= datenum(2016,12,12) || (Time_meas(j) >= datenum(2016,12,16) && Time_meas(j) <= datenum(2016,12,27)) || Time_meas(j) >= datenum(2016,12,31)
        if abs(Time(i) - Time_meas(j)) < 1/24/60

            Time_mod = [ Time_mod ; Time(i)];

            Tp_mod_2 = [Tp_mod_2 ; Tp(i)];
            DO_mod_2 = [DO_mod_2 ; DO(i)];
            pH_mod_2 = [pH_mod_2 ; pH(i)];

            Time_meas_2 = [Time_meas_2; Time_meas(j)];
            Tp_meas_2 = [Tp_meas_2; Tp_meas(j)];
            pH_meas_2 = [pH_meas_2; pH_meas(j)];
            DO_meas_2 = [DO_meas_2; DO_meas(j)];

            j = j + 1;
            i = i+1;
        else
            i = i+1;
        end
    else
        j = j + 1;
    end
    i = i+1;
end


load gong.mat;
sound(y)
%% Screening of data

m = size(Tp_meas_2,1);

Tp_meas_3 = [];
Tp_mod_3 = [];

DO_meas_3 = [];
DO_mod_3 = [];

pH_meas_3 = [];
pH_mod_3 = [];

for i = 1:m
    if Tp_meas_2(i) > 20
        Tp_meas_3 = [Tp_meas_3 ; Tp_meas_2(i)];
        Tp_mod_3 = [Tp_mod_3 ; Tp_mod_2(i)];
    end
end
% 
% figure(300), clf, 
% plot(Tp_meas_3 + 273.15,Tp_mod_3,'xk')
% 
% figure(301), clf, 
% plot(pH_meas_3,pH_mod_3,'xk')
% 
% 
% figure(302), clf, 
% plot(DO_meas_3,DO_mod_3,'xk')
%% Reg Lin

% R_Tp = corrcoef(Tp_meas_2 + 273.15, Tp_mod_2)
% R_pH = corrcoef(pH_meas_2,pH_mod_2)
% R_DO = corrcoef(DO_meas_2,DO_mod_2)
% 
% lm = fitlm(Tp_mod_3,Tp_meas_3 + 273.15,'linear','Intercept',false)
% lm = fitlm(pH_mod_3,pH_meas_3,'linear','Intercept',false)
% lm = fitlm(DO_mod_3,DO_meas_3,'linear','Intercept',false)

lm_Tp = fitlm(Tp_mod_2,Tp_meas_2 + 273.15,'linear','Intercept',false)
lm_Tp.Rsquared.Ordinary

lm_pH = fitlm(pH_mod_2,pH_meas_2,'linear','Intercept',false)
lm_pH.Rsquared.Ordinary

lm_DO = fitlm(DO_mod_2,DO_meas_2,'linear','Intercept',false)
lm_DO.Rsquared.Ordinary

load gong.mat;
sound(y)
%% Figures validation

% break

figure(200), clf, 
plot(Tp_meas_2 ,Tp_mod_2 - 273.15,'xk')

figure(201), clf, 
plot(pH_meas_2,pH_mod_2,'xk')


figure(202), clf, 
plot(DO_meas_2,DO_mod_2,'xk')

figure(210), clf, hold on
    subplot(1,2,2)
    hist(Tp_meas_2 + 273.15 - Tp_mod_2, 1000);
    h = findobj(gca,'Type','patch');
    h.FaceColor = 'k';
    h.EdgeColor = 'k';
    ax = gca; ax.FontSize = 14; 
    xlabel('Pond temperature modelling residuals','FontSize',18), ylabel('Frequency','FontSize',18)

    subplot(1,2,1), hold on
    plot(Tp_mod_2 - 273.15, Tp_meas_2 ,'xk')
    ax = gca; ax.FontSize = 14; 
     xlabel('Predicted pond temperature (°C)','FontSize',18), ylabel('Measured pond temperature (°C)','FontSize',18)
     
    X = 0:1:30; Y = 0:1:30;
    plot(X,Y,'--r','LineWidth',3)
     
figure(211), clf, 
hist(pH_meas_2 - pH_mod_2,1000)
ax = gca; ax.FontSize = 14; 
xlabel('Pond pH modelling residuals','FontSize',18), ylabel('Frequency','FontSize',18)

figure(212), clf, 
hist(DO_meas_2 - DO_mod_2,1000);
ax = gca; ax.FontSize = 14; 
xlabel('Pond DO concentration modelling residuals','FontSize',18), ylabel('Frequency','FontSize',18)

%% 


figure(310), clf, hold on
    subplot(1,2,2)
    hist(Tp_meas_3 + 273.15 - Tp_mod_3, 1000);
    h = findobj(gca,'Type','patch');
    h.FaceColor = 'k';
    h.EdgeColor = 'k';
    ax = gca; ax.FontSize = 14; 
    xlabel('Pond temperature modelling residuals','FontSize',18), ylabel('Frequency','FontSize',18)

    subplot(1,2,1), hold on
    plot(Tp_mod_3 - 273.15, Tp_meas_3 ,'xk')
    ax = gca; ax.FontSize = 14; xlim([min(Tp_mod_3 - 273.15) 30])
     xlabel('Predicted pond temperature (°C)','FontSize',18), ylabel('Measured pond temperature (°C)','FontSize',18)
     
    X = 20:1:30;
    Y = 20:1:30;
    plot(X,Y,'--r','LineWidth',3)
     
    lm_Tp_3 = fitlm(Tp_mod_3,Tp_meas_3 + 273.15,'linear','Intercept',false)
    lm_Tp_3.Rsquared.Ordinary

load gong.mat;
sound(y)