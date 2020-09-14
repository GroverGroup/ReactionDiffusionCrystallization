% this function finds du/dt for all carrier nodes and bulk concentrations and crystal moments
% calls for both reaction kin and crystal kin functions which are needed 
% all units are similar to those in the main code 
               
function dudt = rate_calc(u,drug,carrier,pH_control,pH_set,num,R,v_reactor,e0,N_bead)
Ns = num(1);
N = num(3);
c = u(1 : length(R)*N*Ns+Ns);                                           % all concentrations      
m = u(length(R)*N*Ns+Ns + 1 : end);                                     % all moments
bulks = c(length(R) * (Ns*N) +1 : length(R) * (Ns*N) +Ns);              % bulk values for each species 

if isequal(drug,'CEX')
    mw_cry = 365.4; 
    rho_crys = 1.4e-12;  % crystal density gr/micron^3
elseif isequal(drug,'AMP')
    mw_cry = 403.5; 
    rho_crys = 1.5e-12;  % crystal density gr/micron^3 
end

%% Simulating the control of the bulk pH with adding Na+
Na = bulks(5);
Cl = bulks(6);
if isequal(pH_control,'y')
    Na = titr(bulks,drug,pH_set);   
end 
bulks(5) = Na;
for zz = 1:length(R)                 % updating bulk and interface nodes 
        c(zz*N*Ns - N) = Na;         % simulating adding Na+ to the BULK (neglect volume change, update volume if required)
end                                  % do that for each bead size at the interface node  
                             
%%
%---------------------------------
dcdt = zeros(length(R) * N*Ns + Ns , 1);
dcdt_bulk = zeros(Ns , 1);
for z = 0:length(R) - 1              % now each bead size reacts and interacts with the bulk (internal rxn - diffusion)
    out = RxnDiffusion(drug,carrier,num,c,v_reactor,e0(z+1),N_bead(z+1),R(z+1),length(R),z);
    dcdt( z * (N*Ns) +1 : (z+1) * (N*Ns) , 1) = out(1:N*Ns);
    dcdt_bulk = dcdt_bulk + out(N*Ns + 1 : end );  
end

% dcdt_bulk and interface nodes
dcdt(length(R) * (Ns*N) +1 : end) = dcdt_bulk;
% updating the interface nodes for all bead sizes , should be same as the bulk (equilibrium) 
dcdt(N : N*Ns : end) = dcdt_bulk(1);   % sub
dcdt(2*N : N*Ns : end) = dcdt_bulk(2); % nuc
dcdt(3*N : N*Ns : end) = dcdt_bulk(3); % api
dcdt(4*N : N*Ns : end) = dcdt_bulk(4); % b 
dcdt(5*N : N*Ns : end) = dcdt_bulk(5); % Na
dcdt(6*N : N*Ns : end) = dcdt_bulk(6); % Cl

%---------------------------------
%% calculating the crystallization rate
pH = pH_help([bulks(1:4) ; Na - Cl],drug);
%---------------------------------
% calculation of nucleation and growth rate 
[B, G] = Crystal(drug,bulks,m,num,pH);
%---------------------------------

% only for method of lines
kv = 0.03;                   % volume shape factor of the crystals  - needle 0.03
if (G == 0) 
    n(1) = 0;
else
    n(1) = B/G;
end

dmdt = zeros(num(2),1);
% Calculate the change in moments - method of meoments 
if num(4) == 0
    dmdt(1) = B;
    dmdt(2:5) = G .* (1:4)' .* m(1:4);     
else
    dmdt = [];
end

%% concentration change due to crystallization 
% Calculate the change in solute concentrations (mol/L/min)
% crystallization only for the API species
% note that the decoy node at the interface (i = N) should follow the bulk phase changes
if num(4)==0 
    dcdt(length(R)*N*Ns + 3) = dcdt(length(R)*N*Ns+3) - 3 * kv * G * rho_crys * m(3) / mw_cry;     % API bulk node
    dcdt(N*3) = dcdt(N*3) - 3 * kv * G * rho_crys * m(3) / mw_cry;                                 % API interface node
else % (only for method of lines)
%     dcdt(length(R)*N*Ns+3) =  dcdt(length(R)*N*Ns+3) - 3 * kv * G * rho_crys * sum((bins.^2).*n')/mw_cry * (L(2)-L(1))/num(4);
%     dcdt(N*3) = dcdt(N*3) - 3 * kv * G * rho_crys * sum((bins.^2).*n')/mw_cry * (L(2)-L(1))/num(4);
end
%---------------------------------
dcdt(3*N : N*Ns : end) = dcdt(N*3); % all interface nodes for API should be updated (required only if more than one bead size exists)
%---------------------------------

% Calculate the change in size histogram (only for method of lines)
if num(4)==0
    dndt = zeros(num(4),1);
else
    % Calculate the change in CSD (distribution of crystals in bins)
    %dndx = dss004(Lmin, Lmax, num(1), n(:,i)); % 4th order MOL
    dndx = dss002_modified(L(1), L(2), num(1), n); % 2nd order MOL 
    dndt = -G .* dndx';    
end  

%% time derivative of the system state vector 
dudt = [dcdt ; dmdt ; dndt];

end