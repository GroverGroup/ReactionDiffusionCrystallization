% By Hossein Salami (2020)
% this code is supported by the main function rate_calc which provides the time derivative of the system
% this code deos not consider the crystallization in the beads just in the bulk phase
% vector u contains length(R)*N*Ns element for each species cons. in the bead nodes + Ns bulk conc.
% and +5 additional elements for API crystal moments.
% || A1 A2 A3 A4 Ai | ..., B1 B2 B3 Bi |.... Bulks 1,2,3,4,Ns | Moments 1,2,3,4,5 
% order is always [sub - nuc - API - byprod - Na - Cl]
% Units are: Enzyme mass [gr] - Bead radius [dm] - Volume [Lit] - Time [Min] - Crystal length [micrometer]
%%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++%%%
%%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++%%%

clc
clear
tic

%% Simulation inputs
% [number of species - crystal moments - number of carrier spatial nodes -  CSD bins = 0 for method of moments]
num = [6, 5, 40, 0];
Ns = num(1);                             % number of species
N = num(3);                              % number of spatial nodes in the bead 
                                         % ----------------------------------                                        
drug = 'CEX';                            % 'CEX' or 'AMP'
carrier = 'agarose';                     % carrier TYPE (NOT SIZE) - 'agarose' or 'immobead' 
pH_control = 'y';                        % control pH? ('y' or 'n')
v_reactor = 250e-3;                      % [Lit] reactor volume
enzyme_mass = 20e-3;                     % total [gr] enzyme  
loading = 50e-3;                         % gr enzyme / gr wet bead 
seed_percent = 0.01;                     % percent of seeding (% wt)
tfinal = 100.0;                           % simulation time in minutes
pH_set = 6.5;                            % bulk pH control / initial
s0 = 300e-3;                             % substrate [M] @t=0
n0 = 200e-3;                             % nucleophile [M] @t=0
                                         % ----------------------------------
n_bin = 1;                               % use if more than one R_bead is to be considered 
sample = 1;                              % use to load commercial beads data (1 for ReliZyme, 2 Immobead, etc.)
if n_bin == 1                            % use if want to just use a single Radius value 
    R = 1e-5 * 80;                       % carrier Radius [dm]
    NR = 1;                              % 1: all beads have similar R  
else                                     % for commercial samples with multiple bead sizes (Rs) in the sample 
    Bead_PSD = R_data_preparation(sample,n_bin);
    R = 1e-5 * Bead_PSD(2,:);            % the size range of beads in the sample  
    NR = 0.01 * Bead_PSD(1,:);           % percent in the sample   
end


%%
if isequal(drug,'CEX')
    kw = 10 ^ -14;
    kas = 10 ^ -7.24;
    kan = 10 ^ -6.26;
    kap = 10 ^ -7.12;
    mw_cry = 365.4; 
    mw_API = 347.39;
    rho_crys = 1.4e-12;  % crystal density gr/micron^3
elseif isequal(drug,'AMP')
    kw = 10 ^ -14;
    kas = 10 ^ -7.24;
    kan = 10 ^ -4.83;
    kap = 10 ^ -7.31;
    mw_cry = 403.5; 
    mw_API = 349.5;
    rho_crys = 1.5e-12;  % crystal density gr/micron^3 
end

%% local enzyme concentration calculations
%---------------------------------
enzyme_mol = enzyme_mass / 86500;            % total enzyme mol (using PGA Mw)
beads_mass = enzyme_mass / loading;          % total gr wet bead required  
rho_bead = 1.1e3;                            % gr / dm3 [Lit];
V_bead = (4/3) * pi .* R.^3;                 % single bead volume 
beads_mass = beads_mass .* NR;               % if bead sample has a distribution 
N_bead = beads_mass ./ (V_bead .* rho_bead); % total number of beads (for each bead size if more than one bead radius)

%---------------------------------
% assumption: enzyme loading on each carrier is surface area depedent 
A = enzyme_mass / (4*pi*rho_bead/3) / sum(N_bead .* R.^2);
e0_mass = A ./ R;
e0 = e0_mass * rho_bead / 86500                % local enzyme concentraion in each bead 

etotal = enzyme_mol / v_reactor                % total enzyme concentraion in the reactor 
V_ratio = sum(N_bead .* V_bead) ./ v_reactor;  % Vcat / Vrxr
if abs(sum(N_bead .* V_bead .* e0) - enzyme_mol) > 2e-16  
    disp 'WRONG CALCULATIONS!!!!';
end
%--------------------------------- 

%% Initial concentrations
%---------------------------------
%--------------------------------- INITIAL STATE CONC.
%---------------------------------
p0 = 0e-3;
b0 = 0e-3;
cl0 = s0;                                           % Cl- ion exists as the counter-ion of PGME 
na0 = titr([s0 n0 p0 b0 0 cl0],drug,pH_set);        % Na+ needed to adjust the initial pH
if seed_percent ~= 0                                % required if the solution is seeded 
    out = Csat(drug,[s0 n0 0 b0 na0 cl0],pH_set);
    p0 = out(3);
end
na0 = titr([s0 n0 p0 b0 0 cl0],drug,pH_set);   
 
pH0 = pH_set;
H0 = 10 ^ -pH0;

c0 = zeros( n_bin * (Ns*N) + Ns,1); % vector of initial concentrations (spatial nodes + bulks) 
c0(n_bin * (Ns*N) +1) = s0;         % sub bulk - mol/lit or mol/dm3
c0(n_bin * (Ns*N) +2) = n0;         % nuc bulk
c0(n_bin * (Ns*N) +3) = p0;         % API bulk
c0(n_bin * (Ns*N) +4) = b0;         % byprod bulk
c0(n_bin * (Ns*N) +5) = na0;        % Na bulk
c0(n_bin * (Ns*N) +6) = cl0;        % Cl bulk

% inital concentrations for all internal and interface nodes in all beads
for z = 0 : n_bin-1
c0(z*(N*Ns) + N) = c0(n_bin * (Ns*N) +1);         % sub boundary node (mol/unit volume of bead fluid) - mol/lit
c0(z*(N*Ns)+ 0*N + 1 : z*(N*Ns) + 1*N) = c0(z*(N*Ns) + N);

c0(z*(N*Ns) + 2*N) = c0(n_bin * (Ns*N) +2);       % nuc boundary node
c0(z*(N*Ns)+ 1*N + 1 : z*(N*Ns) + 2*N) = c0(z*(N*Ns) + 2*N);

c0(z*(N*Ns) + 3*N) = c0(n_bin * (Ns*N) +3);       % API boundary node 
c0(z*(N*Ns)+ 2*N + 1 : z*(N*Ns) + 3*N) = c0(z*(N*Ns) + 3*N);

c0(z*(N*Ns) + 4*N) = c0(n_bin * (Ns*N) +4);       % byprod boundary node
c0(z*(N*Ns)+ 3*N + 1 : z*(N*Ns) + 4*N) = c0(z*(N*Ns) + 4*N);

c0(z*(N*Ns) + 5*N) = c0(n_bin * (Ns*N) +5);       % Na boundary node 
c0(z*(N*Ns)+ 4*N + 1 : z*(N*Ns) + 5*N) = c0(z*(N*Ns) + 5*N);

c0(z*(N*Ns) + 6*N) = c0(n_bin * (Ns*N) +6);       % Cl boundary node
c0(z*(N*Ns)+ 5*N + 1 : z*(N*Ns) + 6*N) = c0(z*(N*Ns) + 6*N);
end        

%% initial seeding if any
kv = 0.03;                                % needle crystals shape factor 
initial_seed_mass = (n0) * v_reactor * mw_cry * seed_percent;
load L1.mat                               % example seed crystals CSD
L1 = L1 / 20;                             % simulate ground seeds
[x,z] = hist(L1,40);
L_avg = mean(L1);                         % mean length in micron 
v_rod = kv .* z.^3 ;                      % needle volume for each bin of seeds 
unit_mass = x .* v_rod .* rho_crys;       % needles mass for each bin of seeds 
number_ratio = initial_seed_mass / (sum(unit_mass));  % 
x = x * number_ratio;                   % total number of seed crystals in each bin 
mm0(1,1) = sum(x .* z.^0) / v_reactor;  % 0 moment (#/Lit)
mm0(2,1) = sum(x .* z.^1) / v_reactor;  % 1 moment (micron/Lit)
mm0(3,1) = sum(x .* z.^2) / v_reactor;  % 2 moment (micron^2/Lit)
mm0(4,1) = sum(x .* z.^3) / v_reactor;  % 3 moment ...
mm0(5,1) = sum(x .* z.^4) / v_reactor;  % 4 moment ...

%%
nn0 = zeros(num(4),1);              % no initial crystal bins (use only for method of lines)
u0 = [c0;mm0;nn0];                  % Initial Reactor State Vector 

%%
% uf = []                              % for simulating continuous systems
% q = []                               % input and output flowrates
%%
%--------------------------------- *************
%--------------------------------- Running the dynamic
%--------------------------------- *************
% options = odeset('RelTol',1e-3,'AbsTol',1e-4,'MaxStep',0.05);
options = odeset('MaxStep',0.05);
tspan = [0 tfinal];
[t,u] = ode15s(@(t,u) rate_calc(u,drug,carrier,pH_control,pH_set,num,R,v_reactor,e0,N_bead), tspan, u0, options);

%%
%--------------------------------- 
%---------------------------------  *************
%--------------------------------- Results
%---------------------------------  *************
%--------------------------------- 

L = length(u(:,1));                        % last time point index 
c = u(:,1:length(R)*N*Ns + Ns);            % concentration part of the solution  
moment = u(:,length(R)*N*Ns + Ns+1 : end); % moment part of the solution  
tref = 1;                                  % to normalize the time unit 

subs = c(:,length(R)*(Ns*N)+1);            % bullk values over time (Na needs to be updated for pH control simulation)
nucs = c(:,length(R)*(Ns*N)+2);
apis = c(:,length(R)*(Ns*N)+3);
b = c(:,length(R)*(Ns*N)+4);
na = c(:,length(R)*(Ns*N)+5);
cl = c(:,length(R)*(Ns*N)+6);


%% Batch Process Attributes 
%---------------------------------
[m,ind] = min(nucs);           % max nucleophile conversion point 
Select = apis(ind) / b(ind);
Conver = 100 * (1 - nucs(ind) / nucs(1));
Prod = apis(ind) / (etotal * t(ind) * 60);
disp 'sel,con,prod' , [Select, Conver, Prod]

%% analyzing the results 
% ------
% pH and H of bulk in time (first need to update the Na values for pH control)
pHs = zeros(1,L);
pH0 = pH_set;
for i = 1:L  
    bulks = [subs(i) nucs(i) apis(i) b(i) na(i) cl(i)];
    Na = bulks(5);
    Cl = bulks(6);
    if isequal(pH_control,'y')
        Na = titr(bulks,drug,pH_set);   
    end
na(i) = Na;                      % updating bulk Na
c(i,N*(Ns-1):N*Ns:end) = Na;     % updating all interface nodes Na
pHs(i) = pH_help([subs(i) nucs(i) apis(i) b(i) Na-Cl],drug);      % re-calculating the bulk pH to double check 
end
h = 10 .^ (-pHs);
h = h';
oh = kw ./ h;

% ------
% separated charged states in the bulk in time (if intrested) (A_c -> A + H+)
s = subs ./ (1 + h./kas);
n = nucs ./ (1 + h./kan);
p = apis ./ (1 + h./kap);
s_c = s .* h ./ kas;
n_c = n .* h ./ kan;
p_c = p .* h ./ kap;

% ------
% plot catalyst internal profile at a specific timepoint for each bead
tp = ind; % time step of interest
pHs_pore = zeros(length(R),N);
h_pore = zeros(length(R),N);
for w = 1:length(R)  % loop over all bead sizes             
n_ind = w-1;         % bead size fixed
subs_pore(w,:) = c(tp,n_ind*(N*Ns)+ 0*N + 1 : n_ind*(N*Ns) + 1*N); % substrate pore values in tp time point 
nucs_pore(w,:) = c(tp,n_ind*(N*Ns)+ 1*N + 1 : n_ind*(N*Ns) + 2*N); % nuc pore values in tp time point
apis_pore(w,:) = c(tp,n_ind*(N*Ns)+ 2*N + 1 : n_ind*(N*Ns) + 3*N); % API pore values in tp time point
b_pore(w,:) = c(tp,n_ind*(N*Ns)+ 3*N + 1 : n_ind*(N*Ns) + 4*N);    % byprod pore values in tp time point
na_pore(w,:) = c(tp,n_ind*(N*Ns)+ 4*N + 1 : n_ind*(N*Ns) + 5*N);   % Na pore values in tp time point
cl_pore(w,:) = c(tp,n_ind*(N*Ns)+ 5*N + 1 : n_ind*(N*Ns) + 6*N);   % Cl pore values in tp time point

    for i = 1:N
        pHs_pore(w,i) = pH_help([subs_pore(w,i) nucs_pore(w,i) apis_pore(w,i) b_pore(w,i) na_pore(w,i)-cl_pore(w,i)],drug);
        h_pore(w,i) = 10 ^ (-pHs_pore(i));
    end
    
oh_pore(w,:) = kw ./ h_pore(w,:);
end


%%  concentration-related plots
%---------------------------------
% checking the total charge for electroneutrality
total_charge = n(end) + p(end) + oh(end) + cl(end) - na(end) - h(end) - s_c(end);

% ------------------------
% plots change in the bulk in time
figure(1)
% title('Species bulk profile')
hold on
ts = t/tref;
sty = '-';
plot(ts,subs,sty,'LineWidth',2)
plot(ts,nucs,sty,'LineWidth',2)
plot(ts,apis,sty,'LineWidth',2)
plot(ts,b,sty,'LineWidth',2)
plot(ts,na,'--','LineWidth',2)      % this can be used to calculate how much Na needs to be added by the controller 
legend('PGME','7ADCA','Cex','PG','fontsize',12)
xlabel('Time, min','fontsize',12)
ylabel('Conc., M','fontsize',12)
set(gca,'FontSize',12)
grid on

% ------------------------
% plots distrubutions in pore in a decided timestep above (tp)
figure(2)
title('Distribution inside the bead')
hold on 
dr = R(n_ind+1) / (N-1);
rs = 0:dr:R(n_ind+1);
rs2 = rs * 1e5;
plot(rs2,subs_pore,'-^')
plot(rs2,nucs_pore,'-^')
plot(rs2,apis_pore,'-^')
plot(rs2,b_pore,'-^')
plot(rs2,na_pore,'y-s')
plot(rs2,cl_pore,'k-s')
% plot(rs2,pHs_pore,'-o')
xlabel('r, {\mu}m','fontsize',12)
ylabel('Conc., M','fontsize',12)
set(gca,'FontSize',12)
legend('Sub','Nuc','API','Byprd','Na','Cl')
grid on

%% crystal-related plots
% ------------------------ moments
% plot first moment / number of crystals in time 
figure(3)
hold on
subplot(131)           
plot(ts,moment(:,1))    
xlabel('Time, min','fontsize',12)
ylabel('{\mu}_0, #/lit','fontsize',12)
grid on

% net crystal mass produced in a batch
subplot(132)           
crystal_mass_0 = p0 * v_reactor * mw_cry;    % required crystal to saturate the initial solution for seeding 
crystal_mass = moment(:,4) .* v_reactor * kv * rho_crys -  crystal_mass_0 - initial_seed_mass; % net crystal mass produced 
plot(ts,crystal_mass)
xlabel('Time, min','fontsize',12)
ylabel('Crystal mass, gr','fontsize',12)
set(gca,'FontSize',12)
grid on

% estimation of mean crystal size (micron)
subplot(133)           
plot(ts,moment(:,2) ./ moment(:,1))
xlabel('Time, min','fontsize',12)
ylabel('Length, micron','fontsize',12)
grid on

% ------------------------ 
% plot the bulk supersaturation during the process
for zz = 1:L
    Csature = Csat(drug,c(length(R)*N*Ns+1:length(R)*N*Ns+Ns),pHs(zz));
    ss(zz) = apis(zz) / Csature(3);
end
figure(4)
plot(ts,ss)
xlabel('Time, min','fontsize',12)
ylabel('Bulk supersaturation','fontsize',12)


% ------------------------
% visualizing the evolution of internal bead pH or concentration in time for all bead size groups 
%%% --------------------------------------------------------------------
%%%
% figure(5)
% hold on
% ii = 1;
% for k = 1:length(R)
%     rss(k,:) = 1e5 * linspace(0,R(k),N);
% end
% while t(ii) < 300 % until time 
%       tp = ii;
% for w = 1:length(R)
%     n_ind = w-1; % bead of interest 
%     subs_pore(w,:) = c(tp,n_ind*(N*Ns)+ 0*N + 1 : n_ind*(N*Ns) + 1*N);   % pore values in -tp- time point for -w- bead size group  
%     nucs_pore(w,:) = c(tp,n_ind*(N*Ns)+ 1*N + 1 : n_ind*(N*Ns) + 2*N);
%     apis_pore(w,:) = c(tp,n_ind*(N*Ns)+ 2*N + 1 : n_ind*(N*Ns) + 3*N);
%     b_pore(w,:) = c(tp,n_ind*(N*Ns)+ 3*N + 1 : n_ind*(N*Ns) + 4*N);
%     na_pore(w,:) = c(tp,n_ind*(N*Ns)+ 4*N + 1 : n_ind*(N*Ns) + 5*N);
%     cl_pore(w,:) = c(tp,n_ind*(N*Ns)+ 5*N + 1 : n_ind*(N*Ns) + 6*N);
%     
%         for i = 1:N
%             pHs_pore(w,i) = pH_help([subs_pore(w,i) nucs_pore(w,i) apis_pore(w,i) b_pore(w,i) (na_pore(w,i)-cl_pore(w,i))],drug);
%             h_pore(w,i) = 10 ^ (-pHs_pore(i));
%         end  
%         
% oh_pore(w,:) = kw ./ h_pore(w,:);
% 
% plot(rss(w,:),pHs_pore(w,:),'-o','LineWidth',2)
% % plot(rss(w,:),nucs_pore(w,:))
% 
% end
%      
%     ii = ii + 2000;
%     t(ii);
%     pause(0.2)
% 
% end
% xlabel('r, {\mu}m','fontsize',12)
% ylabel('pH')


% ------------------------
% calculating net crystal mass produced until a fixed PG limit
if max(b) > 0.05
    jj = find(b>0.05, 1 );
else 
    jj = L;
end
time = ts(jj)
mass = crystal_mass(jj)


toc
