function [ K2 ] = K_carbonate_2( Tp )

%% This function aims at calculating the equilibrium constant for the acidic reaction of carbon dioxide and bicarbonate HCO3- <=> HCO32-
% This equilibrium constant is valid for activities calculated in mol/L

% INPUT:

% Tp : Temperature of the pond (K)

Tp = Tp - 273.15;

if Tp < 5
    K2  = 2.754 * 10^(-11);
else
    if 5 < Tp < 10
        K2 = (2.754 + ( 3.236 - 2.754)/(10 - 5)*(Tp - 5))*10^(-11);
    else
        if 10 < Tp < 15
            K2 = (3.236 + ( 3.715 - 3.236)/(15 - 10)*(Tp - 10))*10^(-11);
        else
            if 15 < Tp < 20
                K2 = (3.715 + ( 4.169 - 3.715)/(20 - 15)*(Tp - 15))*10^(-11);
            else
                if 20 < Tp < 25
                    K2 = (4.169 + ( 4.477 - 4.169)/(25 - 20)*(Tp - 20))*10^(-11);
                else
                    if 25 < Tp < 30
                        K2 = (4.477 + ( 5.129 - 4.477)/(30 - 25)*(Tp - 25))*10^(-11);
                    else
                        if 30 < Tp < 40
                            K2 = (5.129 + ( 6.026 - 5.129)/(40 - 30)*(Tp - 30))*10^(-11);
                        else
                            if 40 < Tp 
                                K2 = (6.026)*10^(-11);
                            end
                        end
                    end
                end
            end
        end
    end
end
                        
end

