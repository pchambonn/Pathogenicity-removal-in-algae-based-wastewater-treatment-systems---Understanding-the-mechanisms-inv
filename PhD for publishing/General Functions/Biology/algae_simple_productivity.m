function [ P_algae_gross, P_algae_resp, P_algae_net ] = algae_simple_productivity( Hs_dir, Hs_dif, Tp, sigma, X_algae, d, S, TC, K_C, PE, HV, eps_eff,DO,K_DO_algae,lambda_algae_light,lambda_algae_dark,teta_kd)
% This function aims at calculating the productivity of algae based on a
% simple photosynthetic efficiency coefficient


% INPUTS:

% Hs_dir : direct solar radiation or global if assuming it is only direct (W/m2)
% Hs_dif : diffuse solar radiation (W/m2)
% Tp : pond temperature (K)
% sigma : optical density of the algal medium (m-1)
% X_algae : Algae concentration (kg DW/m3)
% d : pond depth (m)
% S : pond surface (m2)
% TC : Total carbon concentration (g C/m3)
% K_C : affinity of algae on carbon (on Total Carbon mg-C/L, 0.00432 from Solimenio et al. ; 0.035 from Bai et al.)
% PE : Photosynthetic effiency (-)
% HV : algae heat value (J/kg DW)
% eps_eff : safety factors accoutning for extra energy required by algae in non optimal growth conditions (cf Bechet 2013);
% DO : dissolved oxygen concentration (mg O2/l)
% K_DO_algae : affinity of algae on DO during respiration (mg/L) (0.02 from Solimeno et al.)
% lambda_algae_light : rate of respiration of algae in light conditions(kg TSS decayed/kg TSS produced/s)
% lambda_algae_dark : rate of respiration of algae in dark conditions(kg TSS decayed/kg TSS produced/s)
% teta_kd: correction for temperature

% OUTPUTS

% P_algae_gross : algae gross productivity (kg DW/s)
% P_algae_resp : algae mass loss through respiration (kg DW/s)
% P_algae_net : algae net mass productivity (kg DW/s)

    I0_dir = Hs_dir*0.47; % We need to convert the intensity values in PAR from JQ data
    I0_dif = Hs_dif*0.47;
    
    E_used = (I0_dir + I0_dif)*(1 - exp(-sigma*d))*S ; % Power used by algae for photosynthesis during the time step (J/s)
    % The 0.85 factor accounts ofr the fact that algae are in VSS while the
    % attenuation was determined based on TSS (Bechet et al. 2013 for the
    % 15% ash content).
    if E_used > 0.00001
        P_algae_gross = (TC/(K_C + TC))*PE*E_used/(HV*eps_eff); % Algal net production (kg/s)
    
        kd = lambda_algae_light*teta_kd^(Tp - 273.15 - 20);
        P_algae_resp = (DO/(K_DO_algae + DO))*kd*X_algae*d*S;
        P_algae_net = P_algae_gross - P_algae_resp;
    else
        P_algae_gross = 0;
        kd = lambda_algae_dark*teta_kd^(Tp - 273.15 - 20);
        P_algae_resp = (DO/(K_DO_algae + DO))*kd*X_algae*d*S;
        P_algae_net = P_algae_gross - P_algae_resp;
    end
    
end

