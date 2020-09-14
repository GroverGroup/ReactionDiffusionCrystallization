% this function supports the rate_calc and provides the time derivative of concentrations due to reaction 

function dcdt_calc = RxnDiffusion(drug,carrier,num,c,v_reactor,e0,N_bead,R,n,k)
%%
if isequal(drug,'CEX')
    kw = 10 ^ -14;
    kas = 10 ^ -7.24;
    kan = 10 ^ -6.26;
    kap = 10 ^ -7.12;
elseif isequal(drug,'AMP')
    kw = 10 ^ -14;
    kas = 10 ^ -7.24;
    kan = 10 ^ -4.83;
    kap = 10 ^ -7.31;
end

% defining the mesh 
Ns = num(1);
N = num(3);
dr = R / (N-1);
rs = 0:dr:R;
rs = rs';

%------------------------
% arbitary enzyme distribution inside the bead (3 examples of conserved/valid profiles)
% values for the first and last node do not matter 
nn = 0;
if nn == 0 
    % constant enzyme conc
    e_local = e0 * ones(N,1);
elseif nn == 2 
    % ^2 enzyme concentration
    A = 15*e0/14/R^2;
    B = 5*e0/14;
    e_local = A*rs.^2 + B;
elseif nn == 4 
    % ^4 enzyme concentration - all enz in the 50% diameter 
    A = 7*e0/3/R^4;
    e_local = A*rs.^4;   
elseif nn == 8 
    % ^8 enzyme concentration - all enz in the 25% diameter 
    A = 11*e0/3/R^8;
    e_local = A*rs.^8;     
end

%------------------------
% physical parameters for each bead type (Agarose gel, or Immobead/ReliZyme)
if isequal(carrier,'agarose')
    De_s = 5.65e-8;     % dm^2/s
    De_n = 5.71e-8;     % dm^2/s
    De_p = 5.09e-8;     % dm^2/s
    De_b = 5.68e-8;     % dm^2/s
    De_na = 10e-8;      % dm^2/s for sodium 
    De_cl = 15e-8;      % dm^2/s for Cl 
elseif isequal(carrier,'immobead')
    De_s = 2.44e-8;     % dm^2/s
    De_n = 2.19e-8;     % dm^2/s
    De_p = 1.8e-8;      % dm^2/s
    De_b = 2.53e-8;     % dm^2/s
    De_na = 10e-8;      % dm^2/s for Na 
    De_cl = 15e-8;      % dm^2/s for Cl 
end

% charge of the charged state of each species
zs = 1;
zn = -1;
zp = -1;
zna = 1;
zcl = -1;

% total biocatalyst surface area 
A_bead = N_bead * (4 * pi * R^2); 

%------------------------
% for easier handling of the species (when we have multipe bead sizes in the sample)
c_cut = c( k * (N*Ns) +1 : (k+1) * (N*Ns) );

subs = c_cut(1:N);        % s subs in internal nodes
nucs = c_cut(N+1:2*N);    % n nuc in internal nodes
apis = c_cut(2*N+1:3*N);  % p API in internal nodes
b = c_cut(3*N+1:4*N);     % b PG in internal nodes
na = c_cut(4*N+1:5*N);    % Na+ in internal nodes  
cl = c_cut(5*N+1:6*N);    % Cl- in internal nodes  

%% calculate local pHs and charged / uncharged states conc.
%------------------------
local_pH = ones(N,1);
for i = 1:N
    local_pH(i) = pH_help([subs(i) nucs(i) apis(i) b(i) na(i)-cl(i)],drug);
end
h = 10 .^ (-local_pH);

%------------------------
% different states of each species (weak acids/bases)  (A_c --> A + H+)  
s = subs ./ (1 + h./kas);
s_c = s .* h ./ kas;
n = nucs ./ (1 + h./kan);
n_c = n .* h ./ kan;
p = apis ./ (1 + h./kap);
p_c = p .* h ./ kap;

%% define Js for each species (charged or uncharged) - fluxes to nodes
%-----------------------
dsdr = zeros(N,1);
dndr = zeros(N,1);
dpdr = zeros(N,1);
dbdr = zeros(N,1);
dnadr = zeros(N,1);
dcldr = zeros(N,1);
ds_cdr = zeros(N,1);
dn_cdr = zeros(N,1);
dp_cdr = zeros(N,1);

dsdrr = zeros(N,1);
dndrr = zeros(N,1);
dpdrr = zeros(N,1);
dbdrr = zeros(N,1);
dnadrr = zeros(N,1);
dcldrr = zeros(N,1);
ds_cdrr = zeros(N,1);
dn_cdrr = zeros(N,1);
dp_cdrr = zeros(N,1);

% from Nernst-Plank equation 
sigma = zeros(N,1);
dsigma = zeros(N,1);

% finite-diff derivatives - 3pt central for middle points
dsdr(2:N-1) = (s(3:N) - s(1:N-2)) / (2*dr);
dsdrr(2:N-1) = (s(3:N) - 2*s(2:N-1) + s(1:N-2)) / dr^2;
dndr(2:N-1) = (n(3:N) - n(1:N-2)) / (2*dr);
dndrr(2:N-1) = (n(3:N) - 2*n(2:N-1) + n(1:N-2)) / dr^2;
dpdr(2:N-1) = (p(3:N) - p(1:N-2)) / (2*dr);
dpdrr(2:N-1) = (p(3:N) - 2*p(2:N-1) + p(1:N-2)) / dr^2;
dbdr(2:N-1) = (b(3:N) - b(1:N-2)) / (2*dr);
dbdrr(2:N-1) = (b(3:N) - 2*b(2:N-1) + b(1:N-2)) / dr^2;
dnadr(2:N-1) = (na(3:N) - na(1:N-2)) / (2*dr);
dnadrr(2:N-1) = (na(3:N) - 2*na(2:N-1) + na(1:N-2)) / dr^2;
dcldr(2:N-1) = (cl(3:N) - cl(1:N-2)) / (2*dr);
dcldrr(2:N-1) = (cl(3:N) - 2*cl(2:N-1) + cl(1:N-2)) / dr^2;
ds_cdr(2:N-1) = (s_c(3:N) - s_c(1:N-2)) / (2*dr);
ds_cdrr(2:N-1) = (s_c(3:N) - 2*s_c(2:N-1) + s_c(1:N-2)) / dr^2;
dp_cdr(2:N-1) = (p_c(3:N) - p_c(1:N-2)) / (2*dr);
dp_cdrr(2:N-1) = (p_c(3:N) - 2*p_c(2:N-1) + p_c(1:N-2)) / dr^2;
dn_cdr(2:N-1) = (n_c(3:N) - n_c(1:N-2)) / (2*dr);
dn_cdrr(2:N-1) = (n_c(3:N) - 2*n_c(2:N-1) + n_c(1:N-2)) / dr^2;


% 3-pt backwards for last points (N, interface with bulk)
dsdr(N) = (3*s(N) + s(N-2) - 4*s(N-1) ) / (2*dr);
dndr(N) = (3*n(N) + n(N-2) - 4*n(N-1)) / (2*dr);
dpdr(N) = (3*p(N) + p(N-2) - 4*p(N-1)) / (2*dr);
dbdr(N) = (3*b(N) + b(N-2) - 4*b(N-1)) / (2*dr);
dnadr(N) = (3*na(N) + na(N-2) - 4*na(N-1)) / (2*dr);
dcldr(N) = (3*cl(N) + cl(N-2) - 4*cl(N-1)) / (2*dr);
ds_cdr(N) = (3*s_c(N) + s_c(N-2) - 4*s_c(N-1)) / (2*dr);
dn_cdr(N) = (3*n_c(N) + n_c(N-2) - 4*n_c(N-1)) / (2*dr);
dp_cdr(N) = (3*p_c(N) + p_c(N-2) - 4*p_c(N-1)) / (2*dr);

dsdrr(N) = (s(N) - 2*s(N-1) + s(N-2)) / dr^2;
dndrr(N) = (n(N) - 2*n(N-1) + n(N-2)) / dr^2;
dpdrr(N) = (p(N) - 2*p(N-1) + p(N-2)) / dr^2;
dbdrr(N) = (b(N) - 2*b(N-1) + b(N-2)) / dr^2;
dnadrr(N) = (na(N) - 2*na(N-1) + na(N-2)) / dr^2;
dcldrr(N) = (cl(N) - 2*cl(N-1) + cl(N-2)) / dr^2;
ds_cdrr(N) = (s_c(N) - 2*s_c(N-1) + s_c(N-2)) / dr^2;
dn_cdrr(N) = (n_c(N) - 2*n_c(N-1) + n_c(N-2)) / dr^2;
dp_cdrr(N) = (p_c(N) - 2*p_c(N-1) + p_c(N-2)) / dr^2;

% for electrical potential effect (Nernst-Plank eqn)
s1 = (zs*De_s*ds_cdr + zn*De_n*dndr + zp*De_p*dpdr + zna*De_na*dnadr + zcl*De_cl*dcldr);
s2 = (1*De_s*s_c + 1*De_n*n + 1*De_p*p + 1*De_na*na + 1*De_cl*cl);
ds1s2 = (zs*De_s*ds_cdrr + zn*De_n*dndrr + zp*De_p*dpdrr + zna*De_na*dnadrr + zcl*De_cl*dcldrr) .* s2;
ds2s1 = (1*De_s*ds_cdr + 1*De_n*dndr + 1*De_p*dpdr +  1*De_na*dnadr + 1*De_cl*dcldr) .* s1;
sigma = s1 ./ s2;
dsigma = ((ds1s2) - (ds2s1)) ./ (s2).^2;

%%
% transport part of the governing equation 
Ms = De_s * (2*dsdr./rs + dsdrr);
Ms_c = De_s * (2*ds_cdr./rs + ds_cdrr) - De_s*zs * (sigma .* 2.*s_c./rs + ds_cdr .* sigma + s_c .* dsigma);
Mn = De_n * (2*dndr./rs + dndrr) - De_n*zn * (sigma .* 2.*n./rs + dndr .* sigma + n .* dsigma);
Mn_c = De_n * (2*dn_cdr./rs + dn_cdrr);
Mp = De_p * (2*dpdr./rs + dpdrr) - De_p*zp * (sigma .* 2.*p./rs + dpdr .* sigma + p .* dsigma); 
Mp_c = De_p * (2*dp_cdr./rs + dp_cdrr);
Mb = De_b * (2*dbdr./rs + dbdrr); 
Mna = De_na * (2*dnadr./rs + dnadrr) - De_na*zna * (sigma .* 2.*na./rs + dnadr .* sigma + na .* dsigma); 
Mcl = De_cl * (2*dcldr./rs + dcldrr) - De_cl*zcl * (sigma .* 2.*cl./rs + dcldr .* sigma + cl .* dsigma); 

%%
% reaction part of the governing equation 
%--------------------------------------------
% Rxn kinetic 
if isequal(drug,'CEX')
    Ks = 0.38;   % [mol/L] 
    Kp = 0.057;  % [mol/L] 
    k2 = 162;    % [1/s]   
    k3 = 44;     % [1/s] 
    k4 = 316;    % [1/s]   
    km4 = 217;   % [1/s] 
    k5 = 6.3;    % [1/s]   
    % Kn is a function of pH
    PreExp = 0.0005;   % 0.0011 for amp and amx, 0.0005 for cex
    exponent = 0.525;  % 0.525 for amp and amx
    Kn = PreExp*exp(local_pH * exponent);
elseif isequal(drug,'AMP')
    Ks = 0.38;   % [mol/L] 
    Kp = 0.095;  % [mol/L] 
    k2 = 162;    % [1/s]   
    k3 = 44;     % [1/s] 
    k4 = 235;    % [1/s]   
    km4 = 217;   % [1/s] 
    k5 = 9;     % [1/s]   
    % Kn is a function of pH
    PreExp = 0.0011 ;   % 0.0011 for amp and amx
    exponent = 0.525;  % 0.525 for amp and amx
    Kn = PreExp.*exp(local_pH.*exponent);    
end    
Ka1 = 10^-7.52;    % acid dissociation constant (most acidic proton)
Ka2 = 10^-8.19;    % acid dissociation constant (least acidic proton)
% common grouped parameters
A = k3*Kn + k4 * nucs + k5 * nucs;
B = k3*Kn + k5 * nucs;  
% available enzyme at that node
e = e_local ./ ( 1 + subs./Ks + apis./Kp + nucs./Kn + (Kn./A) .* (km4*apis./Kp + k2*subs./Ks) .* ...
                (1 + nucs./Kn + h./Ka1 + Ka2./h) );
           
Rp = e .* ( k2*k4*nucs.*subs ./ (Ks*A) - km4*apis.*B ./ (Kp*A) );
Rb = e .* (B./A) .* (k2*subs/Ks + km4*apis/Kp);
Rn = -Rp;
Rs = -Rb-Rp;
%--------------------------------------------
%%% The governing equation (transport+rxn)
%--------------------------------------------
dsdt = Ms + Ms_c + Rs;
dndt = Mn + Mn_c + Rn;
dpdt = Mp + Mp_c + Rp;
dbdt = Mb + Rb;
dnadt = Mna;
dcldt = Mcl;
%%
% BC at r = 0 (center of the bead)
dsdt(1) = dsdt(2);
dndt(1) = dndt(2);
dpdt(1) = dpdt(2);
dbdt(1) = dbdt(2);
dnadt(1) = dnadt(2);
dcldt(1) = dcldt(2);

% BC at r = R - chemical equilibrium in phase boundary and well-mixed
Js = -De_s * (dsdr(N));
Js_c = -De_s * (ds_cdr(N) - zs * s_c(N) * sigma(N));
Jn = -De_n * (dndr(N) - zn * n(N) * sigma(N));
Jn_c = -De_n * (dn_cdr(N));
Jp = -De_p * (dpdr(N) - zp * p(N) * sigma(N));
Jp_c = -De_p * (dp_cdr(N));
Jb = -De_b * (dbdr(N));
Jna = -De_na * (dnadr(N) - zna * na(N) * sigma(N));
Jcl = -De_cl * (dcldr(N) - zcl * cl(N) * sigma(N));

dsdt(N) = (A_bead / v_reactor) * ( Js + Js_c ) ;
dndt(N) = (A_bead / v_reactor) * ( Jn + Jn_c ) ;
dpdt(N) = (A_bead / v_reactor) * ( Jp + Jp_c ) ;
dbdt(N) = (A_bead / v_reactor) * ( Jb );     
dnadt(N)= (A_bead / v_reactor) * ( Jna );     
dcldt(N)= (A_bead / v_reactor) * ( Jcl );     

dcdt_bulks = zeros(Ns,1);
dcdt_bulks(1) = dsdt(N);  
dcdt_bulks(2) = dndt(N); 
dcdt_bulks(3) = dpdt(N);  
dcdt_bulks(4) = dbdt(N); 
dcdt_bulks(5) = dnadt(N);
dcdt_bulks(6) = dcldt(N);


dcdt = [dsdt ; dndt ; dpdt ; dbdt ; dnadt ; dcldt ; dcdt_bulks];

% checking the net charge flow for Nernst-Plank 
% if (-Jn -Jp - Jcl + Js_c + Jna > 1e-5)
%     display 'net charge flux is too large'
% end

% unit would be [mol] / vol [lit] / time [min]
dcdt_calc = 60 * dcdt;

end
