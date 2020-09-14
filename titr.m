% calculates how much Na+ is required for setting the pH at the set value
% this function is necessary for righ initial conditions and also to calculate what is the required titrant for pH control  
% pKa: [PGME NUC NUC API API PG PG] 

function lhs = titr(c,drug,pH)

if isequal(drug,'AMP')
    pKa_std = [7.24 2.59 4.83 3.26 7.31 2.08 9.14];
elseif isequal(drug,'CEX')
    pKa_std = [7.24 2.59 6.26 3.26 7.12 2.08 9.14];
end
pKa_new = pKa_std + 1;
check = pKa_std;

sub = c(1);
nuc = c(2);
pro = c(3);
byp = c(4);
na = c(5);
cl = c(6);

kw = 10 ^ -14;
h = 10 ^ -pH;
oh = kw / h;
%------------------------
while sum(abs(pKa_new - check)) > 0.01 

check = pKa_new;

kas = 10 ^ (-pKa_new(1));
kan = 10 ^ (-pKa_new(3));
kap = 10 ^ (-pKa_new(5));
kab = 10 ^ (-pKa_new(7));
                                  
s = sub ./ (1 + h./kas);    % calculations of charged states based on A_c -> A + H (with Ka2)
s_c = s .* h ./ kas;
n = nuc ./ (1 + h./kan);
n_c = n .* h ./ kan;
p = pro ./ (1 + h./kap);
p_c = p .* h ./ kap;
b = byp ./ (1 + h./kab);
b_c = b .* h ./ kab;

I = 0.5 * (s_c + n + p + b + cl + na + h + oh);
A = 0.509;            % Debye-Huckel constant (zero if above)
z = [1 1 0 1 0 1 0];  % Charge of the acid

pKa_new = pKa_std + A .* (2.*z-1) .* (sqrt(I)./(1+sqrt(I)) - 0.2.*I);

na = cl + n + p + b + oh - s_c - h;
end

lhs = na;
end
