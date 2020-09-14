% this function supports the rate_calc and provides the crystal birth and growth rate (if supersaturated)
% it calls for function Csat which provides the saturation limit and then calculates the supersaturation based on that 
% please use a valid Csat function which returns solubility of your crystallizing agent in SI units.


function [B, G] = Crystal(drug,bulks,m,num,pH)

Ns = num(1);
N = num(3);
bins = num(4);

sat_conc = Csat(drug,bulks,pH);         % API solubility as a function of pH and Ionic strength  
S_anti = bulks(3)/sat_conc(3);          % API supersaturation

%% crystallization kinetics 
if isequal(drug,'CEX')     
    rho_crys = 1.4e-12;  % crystal density gr/micron^3
    Kj = 2.54;           % [#/min.L]
    B0 = 1.79;           % [-]
    Kb = 2.98 * 1e5;     % [/min.L]
    b = 1;               % [-]
    mm = 0.46;           % [-]
    Kg = 6.52;           % [micron/min]
    g = 2;               % [-]
    
elseif isequal(drug,'AMP')  
    rho_crys = 1.5e-12;  % crystal density gr/micron^3 
    kB1_amp = 5.00e10;   % [#/L/min] 
    B0_amp = 1.27;       % [-]
    kB2_amp = 2.2e9;     % [#/L/min]  
    b_amp = 0.60;        % [-]
    sExp_amp = 1.37;     % [-]
    kG_amp = 8.95;       % [um/min]
    gExp_amp = 1.87;     % [-] 
end

%% calculation of the slurry density 
% For needle like crystals 
kv = 0.03;                   % volume shape factor of the crystals
if num(4) == 0 
    volCryst = m(4) * kv;    % total crystal volume [um^3] per solution volume [lit]  
else
    volCryst = sum((bins.^3).*n') * (L(2)-L(1))/num(4) * kv;
end
% suspension density [g/kg] - [g/Lit]
M = rho_crys * volCryst;

%% calculation of nucleation and growth rate
if isequal(drug,'CEX')
    if (S_anti >=1)
        % Primary nucleation
        B1 = Kj * S_anti * exp( -B0 / (log(S_anti))^2 );
        % Secondary nucleation
        if (M >= 0)
            B2 = Kb * (Kg * (S_anti-1)^g)^b * (M^mm);
        else
            B2 = 0;
        end
        % Growth
        G = Kg * (S_anti-1)^g;

    elseif m(1) < 1     % no supersaturation, no crystal exists 
        B1 = 0;
        B2 = 0;
        G = 0;

    else                % no supersaturation, crystals present --> desolution rate can be implemented
        B1 = 0;
        B2 = 0;
        G = -Kg*1e-5; 
        disp 'disolution of crystals!';
    end
end

%%
if isequal(drug,'AMP')
    if (S_anti >=1)
        % Primary nucleation
        B1 = kB1_amp * exp(-B0_amp / (log(S_anti))^2);
        % Secondary nucleation (if statement handles weird ode45 behavior)
        if (M >= 0)
            B2 = kB2_amp * M^b_amp * (S_anti-1)^sExp_amp;
        else
            B2 = 0;
        end

        % Growth
        G = kG_amp * (S_anti-1)^gExp_amp;

    elseif m(1) < 1     % no super saturation, no crystal exists 
        B1 = 0;
        B2 = 0;
        G = 0;

    else                % no super saturation, crystals present --> desolution rate can be implemented
        B1 = 0;
        B2 = 0;
        G = -kG_amp*1e-5;  
        disp 'disolution of crystals!';
    end
end

% time unit is [min] , crystal length [micron] , solution volume [lit]
B = B1 + B2;

end

