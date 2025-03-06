function [r,p] = alloc_parameters(x,r,p,xi)

r.beta         = x(xi.r_beta);
r.beta_mdr     = r.beta*x(xi.rf_beta_mdr);
r.sym          = x(xi.r_sym);  
p.pu           = x(xi.p_pu);
r.cs           = x(xi.r_cs);   
r.cs2          = x(xi.r_cs)*x(xi.rf_cs2);
r.mort_TB      = r.mort_TB0*x(xi.rf_mort_TB);
p.Dx           = x(xi.p_Dx);
r.default      = r.Tx*(1-x(xi.p_Tx_complete))./x(xi.p_Tx_complete);
p.MDR_rec      = x(xi.p_MDRrec2019);
r.MDR_acqu     = x(xi.r_MDR_acqu);
r.selfcure     = r.selfcure0*x(xi.rf_self_cure);
r.progression  = r.progression0*x(xi.rf_progression);
r.reactivation = r.reactivation0*x(xi.rf_reactivation);
r.relapse      = r.relapse0.*x(xi.rf_relapse);
p.imm          = p.imm0.*x(xi.rf_p_imm); 


