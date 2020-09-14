function Csat = Csat(drug,bulk,pH)

sub = bulk(1);
nuc = bulk(2);
pro = bulk(3);
byp = bulk(4);
na = bulk(5);
cl = bulk(6);

if isequal(drug,'CEX')
    pKa1 = [7.24 2.59 3.26 2.08];  % [PGME, 7ADCA, CEX, PG]
    pKa2 = [7.24 6.26 7.12 9.14];  % [PGME, 7ADCA, CEX, PG]

    % ----- adjustment of pKa ionic strength 
    kw = 10 ^ -14;
    h = 10 ^ -pH;
    oh = kw / h;
    kas = 10 ^ (-pKa2(1));
    kan = 10 ^ (-pKa2(2));
    kap = 10 ^ (-pKa2(3));
    kab = 10 ^ (-pKa2(4));
    s = sub ./ (1 + h./kas);
    s_c = s .* h ./ kas;
    n = nuc ./ (1 + h./kan);
    n_c = n .* h ./ kan;
    p = pro ./ (1 + h./kap);
    p_c = p .* h ./ kap;
    b = byp ./ (1 + h./kab);
    b_c = b .* h ./ kab;
    I = 0.5 * (s_c + n + p + b + cl + na + h + oh);
    A = 0.509;      % Debye-Huckel constant 
    z = [1 0 0 0];  % Charge of the acid
    pKa2 = pKa2 + A .* (2.*z-1).*(sqrt(I) ./ (1+sqrt(I))-0.2.*I);

    SI = 0.038;
    Csat = SI .* (1 + (10.^(-pKa2)) ./ (10.^(-pH)) );
end
%%
if isequal(drug,'AMP')
    pKa = [7.24 4.83 7.31 9.14];  % Santana values for [PGME, 6APA, Amp, PG]
    
    % ----- adjustment of pKa ionic strength
    kw = 10 ^ -14;
    h = 10 ^ -pH;
    oh = kw / h;
    kas = 10 ^ (-pKa(1));
    kan = 10 ^ (-pKa(2));
    kap = 10 ^ (-pKa(3));
    kab = 10 ^ (-pKa(4));
    s = sub ./ (1 + h./kas);
    s_c = s .* h ./ kas;
    n = nuc ./ (1 + h./kan);
    n_c = n .* h ./ kan;
    p = pro ./ (1 + h./kap);
    p_c = p .* h ./ kap;
    b = byp ./ (1 + h./kab);
    b_c = b .* h ./ kab;
    I = 0.5 * (s_c + n + p + b + cl + na + h + oh);
    A = 0.509;  % Debye-Huckel constant (zero if above)
    z = [1 0 0 0];  % Charge of the acid
    pKa = pKa + A .* (2.*z-1).*(sqrt(I) ./ (1+sqrt(I)) - 0.2.*I);

    % Intrinsic protonated solubility (low pH)
    sol1 = [10 0.02 0.0205 0.0350];
    %sol1 = [10 0.02 0.0205 0.0315];
    % Intrinsic deprotonated solubility (high pH)
    sol2 = [10 1 0.0768 0.0633];
    %sol2 = [10 .6 0.0768 02];  % For high PG solubility at pH > 10 (enz no good)
    % Fraction in protonated state
    prot = 1./(1+10.^(pH-pKa));
    % Fraction in deprotonated state
    deprot = (10.^(pH-pKa))./(1+10.^(pH-pKa));

    Csat = min(sol1./prot, sol2./deprot);
end

end

