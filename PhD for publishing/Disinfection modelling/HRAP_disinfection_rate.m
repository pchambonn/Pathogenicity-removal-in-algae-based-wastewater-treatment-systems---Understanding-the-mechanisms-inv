function [ k_s, k_pH,k_20,teta,k_DO ] = HRAP_disinfection_rate( option_dis,d,k_s_m,k_pH_m,k_20_m,teta_m,k_DO_m )

%% This function aims at calculating the components of the disinfection first order decay rate for E. coli in a HRAP according to different authors 

% INPUT
% option_dis : choice for disinfection modelling: 
    % 1 = Marais (1974) 
    % 2 = Mayo  (1995) 
    % 3 = Craggs (2003)
    % 4 = manual adjustment
% d: pond depth (m)
% k_s_m: solar dependent manually adjusted coefficient (m2.W-1.d-1)
% k_pH: pH dependent manually adjusted coefficient (d-1)
% k_20_m: temperature dependent manually adjusted coefficient (d-1)
% teta_m: temperature dependent Arrheniusish manually adjusted coefficient (-)
% K_DO_m: Disolved oxygen concentration dependent manually adjusted coefficient (L.mgDO-1.d-1)
    
% OUTPUT
% kd: first order disinfection rate for E. coli die off (d-1)

if option_dis == 1
    k_s = 0;
    k_pH = 0;
    k_20 = 2.6;
    teta = 1.19;
    k_DO = 0;
 else
     if option_dis == 2
        k_s = 0.001177055*d;
        k_pH = 0.0135;
        k_20 = 0;
        teta = 1;
        k_DO = 0;
     else
         if option_dis == 3
            k_s = 0.083;
            k_pH = 0;
            k_20 = 0.02;
            teta = 1;
            k_DO = 0;
         else
            k_s = k_s_m;
            k_pH = k_pH_m;
            k_20 = k_20_m;
            teta = teta_m;
            k_DO = k_DO_m;
         end
     end
 end

end

