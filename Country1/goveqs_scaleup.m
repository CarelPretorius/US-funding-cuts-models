function out = goveqs_scaleup(t, in, M0, M1, times, i, s, p, sel, agg)
                
% scale = min((t-times(1))/(times(2)-times(1)),1); not correct 
scale = max(min((t-times(1))/(times(2)-times(1)),1),0);

Mt = M1; Mt.lin = M0.lin + scale*(M1.lin-M0.lin);
Mt.Dxlin = M0.Dxlin + scale*(M1.Dxlin-M0.Dxlin);
out = goveqs_basis2(t, in, Mt, i, s, p, sel, agg);
