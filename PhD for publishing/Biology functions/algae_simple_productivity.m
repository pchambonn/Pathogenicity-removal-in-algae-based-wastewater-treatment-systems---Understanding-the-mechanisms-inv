function [ P_algae_gross, P_algae_resp, P_algae_net ] = algae_simple_productivity( Hs_dir, Hs_dif, Tp, sigma, X_algae, d, S, TC, K_C, PE, HV, eps_eff,DO,K_DO_algae,lambda_algae_light,lambda_algae_dark,teta_kd)

% This function aims at calculating the productivity of algae based on a
% simple photosynthetic efficiency coefficient


% INPUTS:

% Hs_dir : direct solar radiation or global if assuming it is only direct (W/m2)
% Hs_dif : diffuse solar radiation (W/m2)
% Tp : pond temperature (K)
% sigma : optical density of the algal medium (m-1)
% X_algae : Algae concentration (kg VSS/m3)
% d : pond depth (m)
% S : pond surface (m2)
% TC : Total carbon concentration (g C/m3)
% K_C : affinity of algae on carbon (on Total Carbon mg-C/L, 0.00432 from Solimeno et al. ; 0.035 from Bai et al.)
% PE : Photosynthetic effiency (-)
% HV : algae heat value (J/kg VSS)
% eps_eff : safety factors accoutning for extra energy required by algae in non optimal growth conditions (cf Bechet et al. 2013);
% DO : dissolved oxygen concentration (mg O2/l)
% K_DO_algae : affinity of algae on DO during respiration (mg/L) (0.02 from Solimeno et al.)
% lambda_algae_light : rate of respiration of algae in light conditions(kg VSS decayed/kg VSS produced/s)
% lambda_algae_dark : rate of respiration of algae in dark conditions(kg VSS decayed/kg VSS produced/s)
% teta_kd: decay correction factor for temperature

% OUTPUTS

% P_algae_gross : algae gross productivity (kg VSS/s)
% P_algae_resp : algae mass loss through respiration (kg VSS/s)
% P_algae_net : algae net mass productivity (kg VSS/s)

    I0_dir = Hs_dir*0.47;                                                            % Total sunlight energy is here converted into PAR sunlight energy
    I0_dif = Hs_dif*0.47;
    
    E_used = (I0_dir + I0_dif)*(1 - exp(-sigma*d))*S ;            % Power used by algae for photosynthesis during the time step (J/s)
                                                                                                    
    if E_used > 0.00001                                                             % Condition checking for light conditions
        P_algae_gross = (TC/(K_C + TC))*PE*E_used/(HV*eps_eff);
    
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

