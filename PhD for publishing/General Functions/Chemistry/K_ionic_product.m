function [ K2 ] = K_ionic_product( Tp )

%% This function aims at calculating the equilibrium constant for the acidic reaction of carbon dioxide and bicarbonate HCO3- <=> HCO32-
% This equilibrium constant is valid for activities calculated in mol/L

% INPUT:

% Tp : Temperature of the pond (K)

Tp = Tp - 273.15;

if Tp < 10
    K2  = 2.92 * 10^(-15);
else
    if 10 < Tp < 20
        K2 = (0.292 + ( 0.681 - 0.292)/(20 - 10)*(Tp - 10))*10^(-14);
    else
        if 20 < Tp < 25
            K2 = (0.681 + ( 1.01 - 0.681)/(25 - 20)*(Tp - 20))*10^(-14);
        else
            if 25 < Tp < 30
                K2 = (1.01 + ( 1.48 - 1.01)/(30 - 25)*(Tp - 25))*10^(-14);
            else
                if 30 < Tp < 40
                    K2 = (1.48 + ( 2.92 - 1.48)/(40 - 30)*(Tp - 30))*10^(-14);
                else
                    if 40 < Tp < 50
                        K2 = (2.92 + ( 5.47 - 2.92)/(50 - 40)*(Tp - 40))*10^(-14);
                    else
                        if 50 < Tp 
                            K2 = (5.47)*10^(-14);
                        end
                    end
                end
            end
        end
    end
end
                        
end

