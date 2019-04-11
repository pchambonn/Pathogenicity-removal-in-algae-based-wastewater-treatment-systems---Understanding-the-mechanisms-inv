function [ Tp,  me, B, Qra_p, Qra_s, Qev, Qconv,  Qcond, Qi, Qr ] = pond_temperature_instant(Tp,Ta,Hs_dir,Hs_dif,RH,qr,Dw_a,Nu,T_soil,Ti,qi,s_t,V,S,L,rho_w,Cp_w,eps_w,sigma_stefan,fa,eps_a,Nu_a,R,Mw,Lw,lambda_a,alpha_a,alpha_s,Ts_ref,ks)
% This function aims at calculating the temperature of a pond after a given
% step of time s_t starting from known conditions


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
%     % Test avec profil de température linéaire
%     Qcond(i) = ks*S*(Ts_ref-Tp(i))/ls_ref ;


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

