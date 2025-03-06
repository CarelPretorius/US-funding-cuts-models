function M = make_model(p, r, i, s, gps)

% --- Get the linear rates ------------------------------------------------
m  = zeros(i.nstates);
m2 = zeros(i.nstates);

for istr = 1:length(gps.strains)
    strain = gps.strains{istr};
    
    gi = @(st) i.(st).(strain);
    Lf  = gi('Lf');
    Ls  = gi('Ls');
    Lfp = gi('Lfp');
    Lsp = gi('Lsp');
    Isc = gi('Isc'); 
    I   = gi('I');     % Symptomatic incidence
    E   = gi('E');
    Dxs = gi('Dx');
    Txs = gi('Tx');
    Tx2s= gi('Tx2');
    Rlo = gi('Rlo');
    Rhi = gi('Rhi');
    R   = gi('R');
    
   
    % Group Selectors
   ismdr   = strcmp(strain, 'mdr');
    
    % --- Fast progression and LTBI stabilisation
    source  = Lf;
    destins =           [Isc,            Ls];
    rates   = [r.progression, r.LTBI_stabil];
    m(destins, source) = m(destins, source) + rates';
    
    % --- Reactivation
    source = Ls; destin = Isc; rate = r.reactivation;
    m(destin, source) = m(destin, source) + rate;
    
      % --- Preventive therapy
    source = Ls; destin = Lsp; rate = r.pt;
    m(destin, source) = m(destin, source) + rate;
    
    source = Lf; destin = Lfp; rate = r.pt;
    m(destin, source) = m(destin, source) + rate;

    
    % --- Fast progression and LTBI stabilisation those who are under preventive therapy
    source  = Lfp; 
    destins = [Isc,                       Lsp]; 
    rates   = [r.progression*(1-0.6),    r.LTBI_stabil];
    m(destins, source) = m(destins, source) + rates';

    % --- Reactivation from those who are under preventive therapy
    source = Lsp; destin = Isc; rate = r.reactivation*(1-0.6);
    m(destin, source) = m(destin, source) + rate;

  
     % --- Sub-clinical to clinical Infection
     source = Isc;
     destin = I;
     rate   = r.sym;                                 % rate of symptom development (estimate this value)
     m(destin, source) = m(destin, source) + rate;   
    
    % --- Primary careseeking, including access to public sector care
    source = I; destin = Dxs.pu; rate = r.cs*(1- r.access)*p.pu ;
    m2(destin, source) = m2(destin, source) + rate;
    
    source = I; destin = Dxs.pr; rate = r.cs*(1-p.pu);
    m2(destin, source) = m2(destin, source) + rate;
    
 % --- Secondary careseeking
    source = E; destin = Dxs.pu; rate = r.cs2*(1- r.access)*p.pu;
    m2(destin, source) = m2(destin, source) + rate;
    
    source = E; destin = Dxs.pr; rate = r.cs2*(1-p.pu);
    m2(destin, source) = m2(destin, source) + rate;
    
        % --- Case finding from subclinical TB
        source = Isc; destin = Dxs.pu; rate = r.cs3*p.pu ;
        m2(destin, source) = m2(destin, source) + rate;
        
        source = Isc; destin = Dxs.pr; rate = r.cs3*(1-p.pu);
        m2(destin, source) = m2(destin, source) + rate;


    
    for ip = 1:length(gps.provs)
        prov = gps.provs{ip};
        
        Dx = Dxs.(prov); Tx = Txs.(prov); Tx2 = Tx2s.(prov);
        
        % --- Diagnosis
        %ismdr   = strcmp(strain, 'mdr');
        pFLinit = p.Dx(ip)*p.Tx_init(ip)*(1 - ismdr*p.MDR_rec(ip)); %*p.Tx_init2(ip)
        pSLinit = p.Dx(ip)*p.Tx_init2(ip)*ismdr*p.MDR_rec(ip);
        p_ltfu  = 1-(pFLinit+pSLinit); %p.Dx(ip)*p.Tx_init(ip);
        
        source  = Dx;
        destins =          [Tx,       Tx2,      E];
        rates   = r.Dx*[pFLinit,  pSLinit,  p_ltfu];
        m(destins, source) = m(destins, source) + rates';
        
        % --- FL Treatment
        pFLcure = p.cure(ip)*(1-ismdr);
        pSLtran = p.SL_trans(ip)*ismdr;
        rMDRacq = r.MDR_acqu*(1-ismdr);
        
        source  = Tx;
        destins = [Rlo            Tx2                        E                              Rhi,             i.Tx.mdr.(prov)];
        rates   = [r.Tx*pFLcure,  r.Tx*(1-pFLcure)*pSLtran,  r.Tx*(1-pFLcure)*(1-pSLtran),  r.default(ip),   rMDRacq];
        m(destins, source) = m(destins, source) + rates';        
        
        % --- SL Treatment
        source  = Tx2;
        destins = [Rlo                 E];          %Rhi];
        rates   = [r.Tx2*p.cure2(ip),  r.Tx2*(1-p.cure2(ip))+ r.default2(ip)];
        m(destins, source) = m(destins, source) + rates';
        
    end    
    
       
    % --- Relapse
    sources = [Rlo, Rhi, R];
    destin  = Isc;
    rates   = r.relapse;
    m(destin, sources) = m(destin, sources) + rates;
    
    sources = [Rlo, Rhi];
    destin  = R;
    rates   = 0.5;
    m(destin, sources) = m(destin, sources) + rates;
    
    % --- Self cure
    source = intersect(s.infectious, s.(strain)); 
    destin = Rhi;
    rate   = r.selfcure;
    m(destin, source) = m(destin, source) + rate;

end
M.lin = sparse(m - diag(sum(m,1)));
M.Dxlin = sparse(m2 - diag(sum(m2,1)));

%keyboard

% --- Get nonlinear rates of TB transmission
U = i.U; inds = [s.Lf, s.Ls,s.Rlo, s.Rhi, s.R];
         indsp = [s.Lfp, s.Lsp];
for istr = 1:length(gps.strains)
    strain = gps.strains{istr};
    m = zeros(i.nstates);
    Lf = i.Lf.(strain);
    m(Lf,  [U, inds]) = 1;
    
    Lfp = i.Lfp.(strain);
    m(Lfp,  indsp) = 1;
    % Immunity
    m(:,inds) = m(:,inds)*(1-p.imm);
    m(:,indsp) = m(:,indsp)*(1-p.imm);
    % Bring all together
    M.nlin.(strain) = sparse(m - diag(sum(m,1)));
end

% --- Getting force-of-infection for DS and DR-TB
m = zeros(2,i.nstates);
m(1,intersect(s.infectious,s.ds))  = r.beta;
m(2,intersect(s.infectious,s.mdr)) = r.beta_mdr;
m(:,s.Isc) = m(:,s.Isc)*0.75;
M.lambda = sparse(m);

% --- Get the mortality rates
m = zeros(i.nstates,2);
m(:,1) = r.mort;
m(s.infectious_wosc,2) = r.mort_TB; 
M.mortvec = m;