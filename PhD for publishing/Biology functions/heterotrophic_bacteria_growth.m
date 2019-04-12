function [ Q_bacteria_gross, Q_bacteria_decay ] = heterotrophic_bacteria_growth( Tp,mu_20,Ks_20,kd_20,teta_mu,teta_Ks,teta_kd,DO,K_DO,bCOD,X_bacteria,d,S )

% This function aims at calculating heterotrophic bacteria productiivty
% following Metcalf and Eddy methodology

% INPUTS
 
% Tp: pond Temperature (K)
% mu_20 : maximal specific growth at 20°C (g VSS/g VSS/s)
% Ks_20 : affinity on bCOD at 20°C (mg bCOD/L)
% kd_20 : decay rate (g VSS/g VSS/s)
% teta_mu : correction coefficient of max specific growth for temperature
% teta_kd : correction coefficient of decay rate for temperature
% teta_Ks : correction coefficient of affinity on bCOD for temperature
% DO : dissolved oxygen concentration (mg O2/l)
% K_DO : affinity of heterotrophic bacteria on DO (mg/L)
% bCOD : biochemical oxygen demand(mg bCOD/L)
% X_bacteria : heterotrophic bacteria initial concentration (g VSS/m3)
% d : pond depth (m);
% S : pond surface (m2);

% OUTPUT

% Q_bacteria_gross : heterotrophic bacteria gross productivity (g VSS/s)
% Q_bacteria_decay : heterotrophic bacteria loss thorugh decay (g VSS/s)

mu = mu_20*teta_mu^(Tp - 273.15 - 20);
Ks = Ks_20*teta_Ks^(Tp - 273.15 - 20);
kd = kd_20*teta_kd^(Tp - 273.15 - 20);

Q_bacteria_gross = (DO/(K_DO + DO))*(mu*bCOD/(Ks + bCOD))*X_bacteria*d*S; 
Q_bacteria_decay = (DO/(K_DO + DO))*kd*X_bacteria*d*S;

end

