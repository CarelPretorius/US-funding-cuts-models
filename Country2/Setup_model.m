clear all;

gps.strains = {'ds','mdr'};
gps.provs = {'pu','pr'};

states0 = {'U'};                                                          
states1 = {'Lf','Ls','Lfp','Lsp','Isc','I','E','Rlo','Rhi','R'};
states2 = {'Dx','Tx','Tx2'};

[i, s, d, lim] = get_addresses({states0}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states1, gps.strains}, i, s, d, lim);
[i, s, d, lim] = get_addresses({states2, gps.strains, gps.provs}, i, s, d, lim);

d = char(d);

s.infectious      = [s.Isc, s.I, s.E, s.Dx, intersect(s.Tx, s.mdr)];
s.infectious_wosc = [s.I, s.E, s.Dx, intersect(s.Tx, s.mdr)];    % Infectious compartment without subclinical (contributing on TB death)
s.prevalent       = unique([s.infectious, [s.Tx,s.Tx2]]);

% --- Include the auxiliaries ---------------------------------------------
auxnames = {'inc', 'notif', 'mort'};
auxinds  = [   2,       3,      1]; 


for ii = 1:length(auxnames)
    inds = lim + [1:auxinds(ii)];
    i.aux.(auxnames{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;


% --- Make aggregators and selectors --------------------------------------
% Selectors for the incidence
%  Incidence: 1.Total, 2.RR-TB
tmp = zeros(2,i.nstates);
tmp(1,s.Isc) = 1;                          
tmp(2,intersect(s.Isc,s.mdr)) = 1;         
agg.inc = sparse(tmp);

tmp = zeros(i.nstates);
tmp(s.Isc,:) = 1;
sel.inc = tmp - diag(diag(tmp));

tmp = zeros(i.nstates);
tmp(intersect(s.Tx,s.mdr),intersect(s.Tx,s.ds)) = 1;
sel.acqu = sparse(tmp - diag(diag(tmp)));

% Selectors for notifications
tmp = zeros(3,i.nstates);
tmp(1,intersect([s.Tx,s.Tx2],s.pu)) = 1;    % public sector notification
tmp(2,intersect([s.Tx,s.Tx2],s.pr)) = 1;    % private sector notification
tmp(3,s.Tx2) = 1;                           % MDR notification
agg.notif = sparse(tmp);

tmp = zeros(i.nstates);
tmp([s.Tx, s.Tx2], s.Dx) = 1;
sel.notif = sparse(tmp - diag(diag(tmp)));


% --- Define variables and ranges -----------------------------------------
xnames = {'r_beta','rf_beta_mdr','r_sym', 'p_pu','r_cs','rf_cs2', 'rf_mort_TB','p_Dx','p_Tx_complete','p_MDRrec2019', 'r_MDR_acqu','rf_self_cure', 'rf_progression', 'rf_reactivation','rf_relapse',  'rf_p_imm'}; 
xnums  =        [1,           1,     1,       1,     1,       1,            1,      2,             2,              2,            1,            1,                1,                1,            3,          1];


xi = []; lim = 0;
for ii = 1:length(xnames)
    inds = lim + [1:xnums(ii)];
    xi.(xnames{ii}) = inds;
    lim = inds(end);
end
% xi.calib = xi.r_mort_TB;
xi.nx = lim;

bds = [];
bds(xi.r_beta,:)          = [0 40];
bds(xi.rf_beta_mdr,:)     = [0 2];
bds(xi.r_sym,:)           = [0.1 100];
bds(xi.p_pu,:)            = [0 1]; %[0 0.9]; %
bds(xi.r_cs,:)            = [0.1 10]; %[0.05 100];
bds(xi.rf_cs2,:)          = [1 40]; %[1 24];
bds(xi.rf_mort_TB,:)      = [0.01 10];
bds(xi.p_Dx,:)            = [0.75 0.9; 0.1  0.7]; % 0.4 0.7];  
bds(xi.p_Tx_complete,:)   = [0.75 0.95; 0.4 0.8];
bds(xi.p_MDRrec2019,:)    = [0  1; 0  0.05];
bds(xi.r_MDR_acqu,:)      = [0  0.06];
bds(xi.rf_self_cure,:)    = [0.75, 1.25];
bds(xi.rf_progression,:)  = [0.75, 1.25];
bds(xi.rf_reactivation,:) = [0.75, 1.25];
bds(xi.rf_relapse,:)      = repmat([0.75, 1.25],3,1);
bds(xi.rf_p_imm,:)        = [0.75, 1.25]; 

prm.bounds = bds';
% --- Define baseline parameter values ------------------------------------

% Natural history
r.progression0  = 0.0826;
r.LTBI_stabil   = 0.872;
r.reactivation0 = 0.0006;
r.selfcure0     = 1/6;
r.mort_TB0      = 1/6;
r.relapse0      = [0.085 0.14 0.0015]; %[0.032 0.14 0.0015];
r.mort          = 1/66;
p.imm0          = 0.5; %0.8;

% Diagnosis stage
r.Dx           = 52;

p.MDR_rec2015 = [0.083  0]; % 
p.Tx_init2 = [0.88 0];
p.SL_trans = [0.88 0]; 

% Treatment stage
p.Tx_init     = [1 1];
r.Tx          = 2;
r.Tx2         = 0.5;   
p.cure        = [1 1];

p.tsrsl     = [0.48  0.01];                                                     % Treatment success rate in pu and pr
r.default2  = r.Tx2*(1-p.tsrsl)./p.tsrsl;
p.cure2     = [0.5 0];  

% Intervention parameter
r.pt          = 0;
r.access      = 0;
r.cs3         = 0;

% --- Bring them all together
prm.p = p; prm.r = r;
ref.i = i; ref.s = s; ref.d = d; ref.xi = xi;
prm.bounds = bds';

% MDR data source 2015-2021
% https://www.who.int/teams/global-tuberculosis-programme/tb-reports/global-tuberculosis-report-2022/tb-disease-burden/2-3-drug-resistant-tb
% --- Get calibration targets ---------------------------------------------
% Country 2 
data.inc2023 = [161 221  291];     % Annual incidence rate 2023
data.mdr2019 = [1.5  2.97  4.6];   % RR TB in 2019
data.mdr2023 = [1.3  2.9  5];      % RR TB in 2023
data.noti2023 = 177*[0.75 1 1.25]; % Annual notification rate 2023
data.mdriniTX = 1.20*[0.5 1 1.5];  % Started second line treatment
data.mort    = [15 26 38];         % Annual mortality rate 2023
data.sym     = [0.36 0.50 0.80];   % proportion of prevalent TB that has symptoms

f1 = get_distribution_fns(data.inc2023,  'lognorm', 0);
f2 = get_distribution_fns(data.mdr2019,  'lognorm', 0);
f3 = get_distribution_fns(data.mdr2023,  'lognorm', 0);
f4 = get_distribution_fns(data.noti2023, 'lognorm', 0);
f5 = get_distribution_fns(data.mort,     'lognorm', 0);
f6 = get_distribution_fns(data.mdriniTX, 'lognorm', 0);
f7 = get_distribution_fns(data.sym,      'lognorm', 0);

lhd.fn  = @(inc2023, mdr2019, mdr2023, noti2023, mort, mdriniTX, sym) 10*f1(inc2023)+ f2(mdr2019)  + 10*f3(mdr2023)+ f4(noti2023)+ f5(mort) +f6(mdriniTX)+ f7(sym); %   
lhd.sgn = -Inf;

save Model_setup;