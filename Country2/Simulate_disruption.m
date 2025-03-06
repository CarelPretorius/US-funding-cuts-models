clear all; 
 
load calibres1; xsto = xsto3; 
[M,I] = max(outsto3); 
load Model_setup;

ix0 = 2e4; nx = 150; dx = round((size(xsto,1)-ix0)/nx); 
xs = xsto(ix0:dx:end,:,1);
pct_xs= prctile(xs,[2.5,50,97.5]);

cs_red = [0.36, 0.36, 0.36];         % Percentage reduction in care-seeking ratein public sector (3 Scenarios)
disruption_dur = [3/12, 3/12, 3/12]; % Disruption period (year) 
recovery_dur = [20, 1, 3/12];        % Recovery period (year)from disruption 

% % --- Get the Regional notifications
tend1   = 2025;  % Disruption starts   
tend2   = 2035+1;  % End date for simulation

% --- Do the simulations 
opts = odeset('NonNegative',[1:i.nstates],'Refine',64,'AbsTol',1e-10,'RelTol',1e-10);


r0 = r; p0 = p;
inct = [];
noti = [];

mk = round(size(xs,1)/25);
for ii = 1:size(xs,1)
    if mod(ii,mk) == 0; fprintf('%0.5g ',ii/mk); end 
    
    [r,p] = alloc_parameters(xs(ii,:),prm.r,prm.p,xi);
 
    % Get the initial conditions
    [out, aux] = obj(xs(ii,:));
    init = aux.soln(end,1:end);
    
    tref = [2023:1:tend2];
   
     % --- Now simulate without intervention (baseline)
    M0 = make_model(p, r, i, s, gps);
 
    geq = @(t,in) goveqs_basis2(t, in, M0, i, s, p, sel, agg); 
    [t0, soln0] = ode15s(geq, [2023 tend2], init, opts);
     
    soln  = soln0;
    t     = t0;
    soln0 = interp1(t, soln, tref);
    
   % Find the initial condition at the starting point for the interventions
    geq = @(t,in) goveqs_basis2(t, in, M0, i, s, p, sel, agg); 
    [ta, solna] = ode15s(geq, [2023 tend1], init, opts);
    initb = solna(end,:);
    
    % --- Next simulate with disruption: scenario1
    r1 = r; p1 = p; 
    r1.access = cs_red(1);
    Dx_red=(95 - ( (1-cs_red(1))*95 + cs_red(1)*60) ) / 95 ; % sens of Xpert =95%, sens of ssm = 60%
    p1.Dx(1) = p.Dx(1)*(1-Dx_red);
     
    M1 = make_model(p1, r1, i, s, gps);
       
    geq = @(t,in) goveqs_scaleup(t, in, M0, M1, tend1 + [0  disruption_dur(1)], i, s, p, sel, agg);
    [tb, solnb] = ode15s(geq, [tend1 tend1+disruption_dur(1)], initb, opts);

    geq = @(t,in) goveqs_scaleup(t, in, M1, M0, tend1+disruption_dur(1) + [0  recovery_dur(1)], i, s, p, sel, agg);
    [tc, solnc] = ode15s(geq, [tend1+disruption_dur(1) tend2], solnb(end,:), opts);

    soln  = [solna; solnb(2:end,:); solnc(2:end,:)]; 
    t     = [ta;    tb(2:end);  tc(2:end)];
    soln1 = interp1(t, soln, tref); 

   % --- Next simulate with disruption: scenario2
    r1 = r; p1 = p; 
    r1.access = cs_red(2);
    Dx_red=(95 - ( (1-cs_red(2))*95 + cs_red(2)*60) ) / 95 ; % sens of Xpert =95%, sens of ssm = 60%
    p1.Dx(1) = p.Dx(1)*(1-Dx_red);
     
    M2 = make_model(p1, r1, i, s, gps);
       
    geq = @(t,in) goveqs_scaleup(t, in, M0, M2, tend1 + [0  disruption_dur(2)], i, s, p, sel, agg);
    [tb, solnb] = ode15s(geq, [tend1 tend1+disruption_dur(2)], initb, opts);

    geq = @(t,in) goveqs_scaleup(t, in, M2, M0, tend1+disruption_dur(2) + [0  recovery_dur(2)], i, s, p, sel, agg);
    [tc, solnc] = ode15s(geq, [tend1+disruption_dur(2) tend2], solnb(end,:), opts);

    soln  = [solna; solnb(2:end,:); solnc(2:end,:)]; 
    t     = [ta;    tb(2:end);  tc(2:end)];
    soln2 = interp1(t, soln, tref); 

       % --- Next simulate with disruption: scenario3
    r1 = r; p1 = p; 
    r1.access = cs_red(3);
    Dx_red=(95 - ( (1-cs_red(3))*95 + cs_red(3)*60) ) / 95 ; % sens of Xpert =95%, sens of ssm = 60%
    p1.Dx(1) = p.Dx(1)*(1-Dx_red);
     
    M3 = make_model(p1, r1, i, s, gps);
       
    geq = @(t,in) goveqs_scaleup(t, in, M0, M3, tend1 + [0  disruption_dur(3)], i, s, p, sel, agg);
    [tb, solnb] = ode15s(geq, [tend1 tend1+disruption_dur(3)], initb, opts);

    geq = @(t,in) goveqs_scaleup(t, in, M3, M0, tend1+disruption_dur(3) + [0  recovery_dur(3)], i, s, p, sel, agg);
    [tc, solnc] = ode15s(geq, [tend1+disruption_dur(3) tend2], solnb(end,:), opts);

    soln  = [solna; solnb(2:end,:); solnc(2:end,:)]; 
    t     = [ta;    tb(2:end);  tc(2:end)];
    soln3 = interp1(t, soln, tref); 

    % --- Bring them all together (baseline, interventions other than vaccine, interventions including vaccine)
   allmat = cat(3, soln0,soln1,soln2,soln3); 
   dallmat = diff(allmat,[],1);
  
    inct(:,:,ii) = squeeze(dallmat(:,i.aux.inc(1),:))*1e5;   % Total incidence
    mort(:,:,ii) = squeeze(dallmat(:,i.aux.mort,:))*1e5; 
    noti(:,:,ii) = squeeze(dallmat(:,i.aux.notif(1),:))*1e5; % Total notification

    inc_incr1 = (sum(inct(1:8,4,:)) - sum(inct(1:8,1,:)))./sum(inct(1:8,1,:));
    inc_incr2 = (sum(inct(1:8,3,:)) - sum(inct(1:8,1,:)))./sum(inct(1:8,1,:));
    inc_incr3 = (sum(inct(1:8,2,:)) - sum(inct(1:8,1,:)))./sum(inct(1:8,1,:));
    
    mor_incr1 = (sum(mort(1:8,4,:)) - sum(mort(1:8,1,:)))./sum(mort(1:8,1,:));
    mor_incr2 = (sum(mort(1:8,3,:)) - sum(mort(1:8,1,:)))./sum(mort(1:8,1,:));
    mor_incr3 = (sum(mort(1:8,2,:)) - sum(mort(1:8,1,:)))./sum(mort(1:8,1,:));
    
end
fprintf('\n');

     indreased_inc1 = prctile(inc_incr1,[2.5,50,97.5])*100; % minimal
     indreased_inc2 = prctile(inc_incr2,[2.5,50,97.5])*100; % moderate
     indreased_inc3 = prctile(inc_incr3,[2.5,50,97.5])*100; % Worst
     
     indreased_mor1 = prctile(mor_incr1,[2.5,50,97.5])*100; % minimal
     indreased_mor2 = prctile(mor_incr2,[2.5,50,97.5])*100; % moderate
     indreased_mor3 = prctile(mor_incr3,[2.5,50,97.5])*100; % Worst

     % -- To estimate the increase in cumulative incidence and mortality
%     [indreased_inc1;indreased_inc2; indreased_inc3]   % Increase in cumulative incidence
%     [indreased_mor1;indreased_mor2; indreased_mor3]   % Increase in cumulative mortality

    
inc_pct  = permute(prctile(inct,[2.5,50,97.5],3),[2,1,3]);
mrt_pct  = permute(prctile(mort,[2.5,50,97.5],3),[2,1,3]);
noti_pct  = permute(prctile(noti,[2.5,50,97.5],3),[2,1,3]);

inct_adj = inc_pct;
mort_adj = mrt_pct;
noti_adj = noti_pct;

outdat.inct = inct_adj;
outdat.mort = mort_adj;
outdat.noti = noti_adj;
save ('country2','outdat')

cols = linspecer(8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ff=figure;
fs = 14; ms = 24; lw = 1.6; ts = 0.1;
pmat1 = cat(4,inct_adj,mort_adj); % (Yearly)  
pmat = permute(pmat1,[3,2,1,4]);    % By 'permute' we can reorganise the matrix-order as required
tis = {'Annual incidence per 100,000','Annual mortality per 100,000'};
subplot(1,2,1);
    for ii = 1:4     
        plt = pmat(:,:,ii,1); 
        pl1(ii,:) = plot(plt(2,:),'Color',cols(ii,:),'linewidth',lw); hold on;
        jbfill(1:size(plt,2),plt(3,:),plt(1,:),cols(ii,:),'None',1,ts); hold on;
     end
     yl = ylim; yl(1) = 0; ylim(yl);
     ylabel(tis{1},'fontsize',16)
     set(gca,'fontsize',fs,'XTick',1:1:size(plt,2),'XTickLabel',2023:2035);
     xtickangle(60)
     legend(pl1(1:4,:),'Baseline','No recovery in funding','Recovery within 1 year','Recovery within 3 months','location','SouthEast','AutoUpdate','off','fontsize',8);
     xlim([1 8]); 
     
     subplot(1,2,2);
    for ii = 1:4    
        plt = pmat(:,:,ii,2); 
        pl1(ii,:) = plot(plt(2,:),'Color',cols(ii,:),'linewidth',lw); hold on;
        jbfill(1:size(plt,2),plt(3,:),plt(1,:),cols(ii,:),'None',1,ts); hold on;
    end
     yl = ylim; yl(1) = 0; ylim(yl);
     ylabel(tis{2},'fontsize',16)
     xlim([1 8]); 
     set(gca,'fontsize',fs,'XTick',1:1:size(plt,2),'XTickLabel',2023:2030);
     xtickangle(60)
     legend(pl1(1:4,:),'Baseline','No recovery in funding','Recovery within 1 year','Recovery within 3 months','location','SouthEast','AutoUpdate','off','fontsize',8);
 sgtitle('Country2')


