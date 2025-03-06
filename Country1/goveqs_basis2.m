function [out, lam] = goveqs_basis2(t, in, M, i, s, p, sel, agg)
                      
invec = in(1:i.nstates);

% Normalise by populations
lam    = (M.lambda*invec/sum(invec));
allmat = M.lin + M.Dxlin + lam(1)*M.nlin.ds + lam(2)*M.nlin.mdr;

out    = allmat*invec;

% Implement deaths
morts = sum(M.mortvec,2).*invec;
out = out - morts;

% Implement births
births = sum(morts);
out(i.U) = out(i.U)+births;

% Get the auxiliaries
out(i.aux.inc)  = agg.inc*(sel.inc.*allmat)*invec;
out(i.aux.notif) = agg.notif*(sel.notif.*allmat)*invec;
out(i.aux.mort) = sum(M.mortvec(:,2).*invec);


% out(i.aux.incslow) = agg.inc*(sel.incslow.*allmat)*invec;
% out(i.aux.incfast) = agg.inc*(sel.incfast.*allmat)*invec;
