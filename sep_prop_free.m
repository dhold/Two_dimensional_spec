%[rho_gg,rho_eg,rho_ge,rho_ee,rho_fe,rho_ef,rho_ff,rho_fg,rho_gf]
function [rho_out,t_range] = sep_prop_free(t0,tend,np,rho_init,tol)

persistent prop_gg  prop_eg  prop_ge prop_ee prop_fe prop_ef prop_ff prop_fg prop_gf
persistent tmp1 tmp1a tmp1b tmp2 tmp2a tmp2b tmp3 tmp3a tmp3b
if iscell(t0)
    
 [tmp1,tmp1a, tmp1b, tmp2,tmp2a, tmp2b,tmp3,tmp3a,tmp3b] =  deal(t0{:});
 rho_out = []; t_range = []; 

prop_op_full = tend; om_s = np; clear tend
%for scaling coherences   
scale1 = -1i*speye(sum(tmp1a))*om_s; 
scale2 = -1i*speye(sum(tmp2a))*om_s;  
scale3 = -2i*speye(sum(tmp3a))*om_s;
    
prop_op_gg = prop_op_full(tmp1,tmp1); prop_gg = @(t,v) prop_op_gg*v;
prop_op_eg = prop_op_full(tmp1a,tmp1a); prop_eg = @(t,v) (prop_op_eg-scale1)*v;
prop_op_ge = prop_op_full(tmp1b,tmp1b); prop_ge = @(t,v) (prop_op_ge+scale1)*v;
prop_op_ee = prop_op_full(tmp2,tmp2); prop_ee = @(t,v) prop_op_ee*v;
prop_op_fe = prop_op_full(tmp2a,tmp2a); prop_fe = @(t,v) (prop_op_fe-scale2)*v;
prop_op_ef = prop_op_full(tmp2b,tmp2b);  prop_ef = @(t,v) (prop_op_ef+scale2)*v;
prop_op_ff = prop_op_full(tmp3,tmp3);  prop_ff = @(t,v) prop_op_ff*v;
prop_op_fg = prop_op_full(tmp3a,tmp3a);  prop_fg = @(t,v) (prop_op_fg-scale3)*v;
prop_op_gf = prop_op_full(tmp3b,tmp3b);  prop_gf = @(t,v) (prop_op_gf+scale3)*v;

return
end
t_range = linspace(t0,tend,np); 
rho_out = zeros(length(t_range),length(rho_init));

%calculate each component
rho_0 = rho_init(tmp1);  rho_out(:,tmp1) = OD_wrapper(t_range,prop_gg,rho_0,[],'ode45',tol);
rho_0 = rho_init(tmp1a); rho_out(:,tmp1a)= OD_wrapper(t_range,prop_eg,rho_0,[],'ode45',tol);
rho_0 = rho_init(tmp1b); rho_out(:,tmp1b) = OD_wrapper(t_range,prop_ge,rho_0,[],'ode45',tol);
rho_0 = rho_init(tmp2);  rho_out(:,tmp2) = OD_wrapper(t_range,prop_ee,rho_0,[],'ode45',tol);
rho_0 = rho_init(tmp2a); rho_out(:,tmp2a) = OD_wrapper(t_range,prop_fe,rho_0,[],'ode45',tol);
rho_0 = rho_init(tmp2b); rho_out(:,tmp2b) = OD_wrapper(t_range,prop_ef,rho_0,[],'ode45',tol);
rho_0 = rho_init(tmp3);  rho_out(:,tmp3) = OD_wrapper(t_range,prop_ff,rho_0,[],'ode45',tol);
rho_0 = rho_init(tmp3a); rho_out(:,tmp3a) = OD_wrapper(t_range,prop_fg,rho_0,[],'ode45',tol);
rho_0 = rho_init(tmp3b); rho_out(:,tmp3b) = OD_wrapper(t_range,prop_gf,rho_0,[],'ode45',tol);

end
