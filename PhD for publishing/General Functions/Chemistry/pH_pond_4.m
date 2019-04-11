function [ pH ] = pH_pond_4( IC , IP, IN, Sigma, Tp, KP2, KP3, KN )

% This function aims at calculating the pH of a solution depending on a
% previous equilibrium, a new amount of inorganic carbon and a new supply
% of inorganic phosphate. This option is based on the methodology outlined
% in the thesis manuscript.

% INPUTS

% TC : total inorganic carbon (g-C/m3)
% TP : total inorganic phosphorous (g-P/m3)
% IN : total ammoniac N (g-N/m3)
% Sigma : balance of inert ions (mol eq/L)
% Tp : pond temperature (K)
% KP: equilibrium constants of  H2PO4- <=> HPO4(2-) (other phosphate
% equilibria are neglected due to pH between expected between 7 - 11)
% KN: equilibrium constants of  NH4+- <=> NH3

% OUTPUTS

% pH :  calculated pH of the pond

KC1 = K_carbonate_1(Tp);
KC2 = K_carbonate_2(Tp);
Kw = K_ionic_product(Tp);

% Conversion of IP, IN and IC to mol/L
IC = IC/12/1000;
IP = IP/31/1000;
IN = IN/14/1000;

a_pol = Sigma + KC1 + KP2 + KN + IN - IP;
b_pol = Sigma*(KC1 + KP2 + KN) + KC1*KP2 + KC1*KN + KP2*KN + KP2*KP3 + KC1*KC2 + IN*(KC1 + KP2) - Kw - KC1*IC - IP*(KC1 + KN) - 2*IP*KP2;
c_pol = Sigma*( KC1*KP2 + KC1*KN + KP2*KN + KC1*KC2 + KP2*KP3) + KC1*KC2*KP2 + KC1*KC2*KN + KC1*KN*KP2 + KP2*KP3*KN + KP2*KP3*KC1 + IN*(KC1*KP2 + KC1*KC2 + KP2*KP3) - Kw*(KC1 + KP2 + KN) - KC1*IC*(KP2 + KN) - 2*KC1*KC2*IC - IP*(KC1*KC2 + KC1*KN) - 2*IP*KP2*(KC1 + KN) - 3*IP*KP2*KP3;
d_pol = Sigma*(KC1*KC2*KP2 + KC1*KC2*KN + KC1*KN*KP2 + KC1*KP2*KP3 + KN*KP2*KP3) + KC1*KC2*KP2*KN + KC1*KC2*KP2*KP3 + KC1*KN*KP2*KP3 + IN*(KC1*KC2*KP2 + KC1*KP2*KP3) - Kw*(KC1*KP2 + KC1*KN + KP2*KN + KC1*KC2 + KP2*KP3) - KC1*IC*(KN*KP2 + KP2*KP3) - 2*KC1*KC2*IC*(KP2 + KN) - IP*KC1*KC2*KN - 2*IP*KP2*(KC1*KC2 + KC1*KN) - 3*IP*KP2*KP3*(KC1 + KN);
e_pol = KC1*KC2*KP2*KP3*KN + Sigma*(KC1*KC2*KP2*KN + KC1*KC2*KP2*KP3 + KC1*KN*KP2*KP3) - Kw*(KC1*KC2*KP2 + KC1*KC2*KN + KC1*KN*KP2 + KN*KP2*KP3 + KC1*KP2*KP3) + IN*KC1*KC2*KP2*KP3- 2*KC1*KC2*IC*(KN*KP2 + KP2*KP3) - 2*IP*KP2*KC1*KC2*KN - 3*IP*KP2*KP3*(KC1*KC2 + KC1*KN) - KC1*IC*KN*KP2*KP3 ;
f_pol = Sigma*(KC1*KC2*KN*KP2*KP3) - Kw*(KC1*KC2*KP2*KN + KC1*KN*KP2*KP3 + KC1*KC2*KP2*KP3) - 3*IP*KP2*KP3*KC1*KC2*KN - 2*KC1*KC2*IC*KN*KP2*KP3;
g_pol = - Kw*KC1*KC2*KN*KP2*KP3;

polynomial = [g_pol f_pol e_pol d_pol c_pol b_pol a_pol 1];
r_pol = roots(polynomial);

if r_pol(1) > 0 
    nh = 1/real(r_pol(1));
else
    if r_pol(2) > 0 
        nh = 1/real(r_pol(2));
    else
        if r_pol(3) > 0 
            nh = 1/real(r_pol(3));
        else
            if r_pol(4) > 0 
                nh = 1/real(r_pol(4));
            else
                if r_pol(5) > 0;
                    nh = 1/real(r_pol(5));
                else
                    if r_pol(6) > 0;
                        nh = 1/real(r_pol(6));
                    else
                        if r_pol(6) > 0;
                            nh = 1/real(r_pol(7));
                        else
                            nh = 10^(-11);
                        end
                    end
                end
            end
        end
    end
end

pH = -log10(nh);

end



