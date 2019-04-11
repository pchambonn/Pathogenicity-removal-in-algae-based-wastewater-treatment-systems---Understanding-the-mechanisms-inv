function [ K1 ] = K_carbonate_1( Tp )

%% This function aims at calculating the equilibrium constant for the acidic reaction of carbon dioxide and bicarbonate H2CO3* <=> HCO3-
% This equilibrium constant is valid for activities calculated in mol/L

% INPUT:

% Tp : Temperature of the pond (K)

Tp = Tp - 273.15;

if Tp < 5
    K1  = 3.020 * 10^(-7);
else
    if 5 <= Tp < 10
        K1 = (3.020 + ( 3.467 - 3.020)/(10 - 5)*(Tp - 5))*10^(-7);
    else
        if 10 <= Tp < 15
            K1 = (3.467 + ( 3.802 - 3.467)/(15 - 10)*(Tp - 10))*10^(-7);
        else
            if 15 <= Tp < 20
                K1 = (3.802 + ( 4.169 - 3.802)/(20 - 15)*(Tp - 15))*10^(-7);
            else
                if 20 <= Tp < 25
                    K1 = (4.169 + ( 4.467 - 4.169)/(25 - 20)*(Tp - 20))*10^(-7);
                else
                    if 25 <= Tp < 30
                        K1 = (4.467 + ( 4.677 - 4.467)/(30 - 25)*(Tp - 25))*10^(-7);
                    else
                        if 30 <= Tp < 40
                            K1 = (4.677 + ( 5.012 - 4.677)/(40 - 30)*(Tp - 30))*10^(-7);
                        else
                            if 40 < Tp 
                                K1 = (5.012)*10^(-7);
                            end
                        end
                    end
                end
            end
        end
    end
end
                        
end

