
function [out, aux] = get_objective_wRNTCP(x, prm, ref, sel, agg, gps, calfn)
                                         
r = prm.r; p = prm.p; i = ref.i; s = ref.s; xi = ref.xi;

% First, check whether the proposed parameter ranges lie within the
% specified ranges (i.e. implement the uniform priors - can also be done
% after performing the simulations and alongside the likelihood, but this
% way saves time, rejecting out-of-range parameters without needing to do
% the simulations

mat = [prm.bounds(2,1:length(x))-x; x-prm.bounds(1,1:length(x))];
% All values should be positive if x is within the specified ranges

if min(mat(:)) < 0
    
    out = Inf*calfn.sgn;
    % +Inf if using sum of squares, -Inf if using likelihoods
    aux = NaN;
%     keyboard;
else
    
    % Note the use of the subfunction alloc_parameters, to keep things tidy
    % with the allocation of parameters
    [r,p] = alloc_parameters(x,r,p,xi);
      
 % --- Set up the necessary models -----------------------------------------

% Final conditions
p3 = p; r3 = r;
%p3.MDR_rec = p.MDR_rec2019;
M3 = make_model(p3, r3, i, s, gps);

% MDR_rec in 2015
p2 = p; r2 = r;
p2.MDR_rec = p.MDR_rec2015;
M2 = make_model(p2, r2, i, s, gps);

% In absence of RNTCP
p1 = p; r1 = r;
p1.pu = 0; p1.MDR_rec  = [0 0]; p1.cure2 = [0 0]; 
M1 = make_model(p1, r1, i, s, gps);

% In absence of MDR and RNTCP
p0 = p1; r0 = r1;
p0.pu = 0; r0.beta_mdr = 0; r0.MDR_acqu = 0;
M0 = make_model(p0, r0, i, s, gps);

% --- Solve the models ----------------------------------------------------

% Equilibrium model
init = zeros(1,i.nx); seed = 1e-6; init(i.U) = (1-seed); init(i.I.ds) = seed;
geq = @(t,in)goveqs_basis2(t, in, M0, i, s, p0, sel, agg);
[t0, soln0] = ode15s(geq, [0:2e3], init, odeset('NonNegative',[1:i.nstates]));


% Introduce MDR-TB
init = soln0(end,:);
geq = @(t,in)goveqs_basis2(t, in, M1, i, s, p1, sel, agg);
[t1, soln1] = ode15s(geq, [1970 1997], init, odeset('NonNegative',[1:i.nstates]));

% RNTCP scale-up
init = soln1(end,:);
[t2, soln2] = ode15s(@(t,in) goveqs_scaleup(t, in, M1, M2, [1997 2007], i, s, p, sel, agg), [1997 2015], init, odeset('NonNegative',[1:i.nstates]));

% DST scale-up (Assuming from 2015 to 2023)
init = soln2(end,:);
[t3, soln3] = ode15s(@(t,in) goveqs_scaleup(t, in, M2, M3, [2015 2023], i, s, p, sel, agg), [2015 2024], init, odeset('NonNegative',[1:i.nstates]));

%keyboard;

allt   = [t2; t3(2:end)]; % ; t3(2:end); t4(2:end); t5(2:end)];
allsol = [soln2; soln3(2:end,:)]; %; soln3(2:end,:); soln4(2:end,:); soln5(2:end,:)];

 soln = interp1(allt, allsol, 1997:2024);  
%soln  = interp1(t2, soln2, 1997:2020);
dsol   = diff(soln,[],1);

     sfin    = soln(end,:);
%     Nfin    = sum(sfin(1:i.nstates)); fac = 1e5/Nfin;
%     dsol    = diff(soln,1);

  % --- Get the objectives ----------------------------------------------
        incd = dsol(:,i.aux.inc(1))*1e5;
        mdr  = dsol(:,i.aux.inc(2))*1e5;
        mort  = dsol(:,i.aux.mort)*1e5;
        noti   = dsol(:,i.aux.notif(1))*1e5;
        mdriniTX   = dsol(:,i.aux.notif(3))*1e5;
    
       inc2023 = incd(end); %dsol(t2==2019, i.aux.inc(1))*1e5;
       mdr2019 = mdr(end-4); %mdr(end-8);
       mdr2023 = mdr(end); %dsol(t2==2019, i.aux.inc(2))*1e5;
       mort2023= mort(end); %dsol(t2==2019, i.aux.mort)*1e5;
       noti2023= noti(end); %dsol(t2==2019, i.aux.notif(1))*1e5;
       mdriniTX2023= mdriniTX(end); %dsol(t2==2019, i.aux.notif(3))*1e5;
       sym_sim =  sum(sfin(s.I))/sum(sfin(s.prevalent)); 

        
    % Compose the objective function
    out = calfn.fn(inc2023, mdr2019, mdr2023, noti2023,mort2023,mdriniTX2023,sym_sim);
  
  
    % --- Get additional outputs and package ------------------------------
    
    aux.soln  = soln; %[soln2, t2];
    %aux.soln_test = soln_test;
    aux.incd = incd;
    aux.mdr = mdr;
    aux.noti = noti;
    aux.mort = mort;
    aux.mdriniTX = mdriniTX;
    aux.sym      = sym_sim;
  
%     keyboard;
end

end
