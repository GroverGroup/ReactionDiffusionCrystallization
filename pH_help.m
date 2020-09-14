function [pHfunc] = pH_help(c,drug)

%pKa: [PGME NUC NUC API API PG PG] 
if isequal(drug,'AMP')
    pKa = [7.24 2.59 4.83 3.26 7.31 2.08 9.14];
elseif isequal(drug,'CEX')
    pKa = [7.24 2.59 6.26 3.26 7.12 2.08 9.14];
end

pKa_std = pKa;

% Need to iterate to adjust the pKa as the ionic strength changes with pH
% which changes the pKa which changes the pH etc.
sse = 1;
while sse > 0.01

    Ka = 10.^-pKa;
    Kw = 10^-14;

    % AH2+ <--> AH + H+ <--> A- + 2H+
    % Concentration of AH2+ in above equilibrium
    AH2pos = @(H,Ka1,Ka2,TA) (TA*H^2) / (H^2 + H*Ka1 + Ka1*Ka2);

    % Concentration of A- in above equilibrium
    Aneg = @(H,Ka1,Ka2,TA) (TA*Ka1*Ka2) / (H^2 + H*Ka1 + Ka1*Ka2);


    chargeBalance = @(H) c(5) + H - Kw/H + c(1)*H/(H+Ka(1)) +...
        AH2pos(H,Ka(2),Ka(3),c(2)) - Aneg(H,Ka(2),Ka(3),c(2)) +...
        AH2pos(H,Ka(4),Ka(5),c(3)) - Aneg(H,Ka(4),Ka(5),c(3)) +...
        AH2pos(H,Ka(6),Ka(7),c(4)) - Aneg(H,Ka(6),Ka(7),c(4));

    % Newton-Raphson was found to consistantly converge to the same solution
    % across many concentrations of acid or base.  MATLAB's built in solvers
    % will sometimes converge to different minima or not converge to realistic
    % solutions.  Newton-Raphson requires the derivative of the charge balance.

    % Derivative of the concentration of AH2+ in charge balance
    dAH2posdH = @(H,Ka1,Ka2,TA) (TA*H*Ka1*(H+2*Ka2))/(H^2 + H*Ka1 + Ka1*Ka2)^2;
    % Derivative of the concentration of A- in charge balance
    dAnegdH = @(H,Ka1,Ka2,TA) (TA*Ka1*Ka2*(2*H+Ka1))/(H^2 + H*Ka1 + Ka1*Ka2)^2;
    % Derivative of the charge balance with respect to Proton concentration
    dCBdH = @(H) 1 + Kw/H^2 + c(1)*Ka(1)/(H+Ka(1))^2 +...
        dAH2posdH(H,Ka(2),Ka(3),c(2)) + dAnegdH(H,Ka(2),Ka(3),c(2)) +...
        dAH2posdH(H,Ka(4),Ka(5),c(3)) + dAnegdH(H,Ka(4),Ka(5),c(3)) +...
        dAH2posdH(H,Ka(6),Ka(7),c(4)) + dAnegdH(H,Ka(6),Ka(7),c(4));

    % Use Newton-Raphson to find the concentration of protons
    x0 = 10^-14;
    xi = x0 - chargeBalance(x0)/dCBdH(x0);
    xi2 = xi - chargeBalance(xi)/dCBdH(xi);
    while abs(xi2 - xi) > 1e-15
        xi = xi2;
        xi2 = xi2 - chargeBalance(xi2)/dCBdH(xi2);
    end
    HplusNR = xi2;

    % Need to correct the pKas based on ionic strength
    IonicStrength = 0.5*(c(5) + HplusNR + Kw/HplusNR + c(1)*HplusNR/(HplusNR+Ka(1)) +...
        AH2pos(HplusNR,Ka(2),Ka(3),c(2)) + Aneg(HplusNR,Ka(2),Ka(3),c(2)) +...
        AH2pos(HplusNR,Ka(4),Ka(5),c(3)) + Aneg(HplusNR,Ka(4),Ka(5),c(3)) +...
        AH2pos(HplusNR,Ka(6),Ka(7),c(4)) + Aneg(HplusNR,Ka(6),Ka(7),c(4)));
    
    % Dissociations of type HA+ <-> H+ + A do not count because both have a
    % single positively charged species on each side
    A = 0.509;  % Debye-Huckel constant (zero if above)
    z = [1 1 0 1 0 1 0];  % Charge of the acid
    pKa_new = pKa_std + A.*(2.*z-1).*(sqrt(IonicStrength)./...
        (1+sqrt(IonicStrength))-0.2.*IonicStrength);
    sse = sum((pKa_new-pKa).^2);
    pKa = pKa_new;
end

pHfunc = -1*log10(HplusNR);
Is = IonicStrength;
end
