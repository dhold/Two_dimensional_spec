function [rho_eg,rho_door_ee,rho_door_gg]=doorway_fn_HEOM(rho_0,N,Nv,HL,H_trunc_1,...
                        t1_range,tp_max,om_scale,sup_op_ge,sup_op_gg,sup_op_ee,...
                        E1,E2,omega_1,omega_2,tol)
%  om_scale = typical frequency of exciton transfer
%  L_eg is -1i[H,rho_{eg}], plus extra acting on one specific tier of the
%  density matrix set of the HEOM.  But only the bit which mixes this block
%  onto itself obviously.  H_mix_op has the mixing a decay of higher tiers
%   - sum_j Trunc_op_j [V_j,[V_j ,rho]] - n*gamma_n  rho_n

if ~exist('tol','var')
    tol = [1e-7,1e-5]; %tolerances
end 
if isempty(E2) && isempty(omega_2) %pass empty for pump probe
    E2 = conj(E1); omega_2 = NaN; ispumpprobe = true;
    rho_door_ee = zeros(N^2*Nv^2*H_trunc_1,length(omega_1),N,N);
    rho_door_gg = zeros(Nv^2*H_trunc_1,length(omega_1),N,N);
else
    ispumpprobe = false;
    rho_door_ee = zeros(N^2*Nv^2*H_trunc_1,length(omega_1),length(omega_2),N,N);  
    rho_door_gg = zeros(Nv^2*H_trunc_1,length(omega_1),length(omega_2),N,N);
end
    
    rho_eg = zeros(N*Nv^2*H_trunc_1,length(t1_range),N);    
    
for j = 1:N
    
   %this uses the arrow convention of "Advancing HEOM for efficient
   %evalution of ..." 1) is a right arrow 2) is a left arrow
   %V_eg = zeros(N,1);  V_ge = zeros(1,N); V_eg(j) = 1;  V_ge(j) = 1; 
    VV_eg = sparse((j-1)*Nv+1:j*Nv,1:Nv,ones(Nv,1),N*Nv,Nv); 
    VV_ge = sparse(1:Nv,(j-1)*Nv+1:j*Nv,ones(Nv,1),Nv,N*Nv);

   V_eg_right{j} = kron(speye(N*Nv),VV_eg);
   V_eg_right_HEOM{j} = kron(speye(HL),V_eg_right{j});
   V_eg_left{j} = kron(VV_eg.',speye(Nv));     
   V_eg_left_HEOM{j} = kron(speye(HL),V_eg_left{j});
   
   V_ge_right{j} = kron(speye(Nv),VV_ge);
   V_ge_left{j} = kron(VV_ge.',speye(N*Nv));        
   V_eg_right_HEOM{j} = kron(speye(HL),V_ge_right{j}); 
   V_eg_left_HEOM{j} = kron(speye(HL),V_ge_left{j});    
    

    om_ex = om_scale(j);%typical frequency of exciton transition
    
    sup_op_ge_sc = sup_op_ge + 1i*eye(N*Nv^2*HL)*om_ex; %scale with this freq
    rho_eg_prop = @(t,v) sup_op_ge_sc*v;   %equation to solve
    
rho_init = VV_eg*rho_0(1:Nv,1:Nv); %jth exciton transition
%this is not trivial as the HEOM code works in the site basis 
out_trunc = numel(rho_init)*HL;  %no truncation yet
%reshape to column vector and pad with zeros to full size
if HL>1
rho_init = [reshape(rho_init,numel(rho_init),1);repmat(zeros(size(rho_init)),[HL-1,1])]; 
end
%solve the resulting ODE with the solver
rho_eg_tmp = (OD_wrapper(t1_range,rho_eg_prop,rho_init,out_trunc,'ode45',tol));
rho_eg(:,:,j) =  rho_eg_tmp;

end

%% Next calculate the doorway function

    rho_gg_prop = @(t,v) sup_op_gg*v;  
    rho_ee_prop = @(t,v) sup_op_ee*v; 
    
for j1= 1:N %first exciton transition
    om_ex = om_scale(j1);
    for j2 = 1:N %second exciton transition
          
      %  if method == 1
        %integrate over every time step
        for lp = 1:length(t1_range)
            t1 = t1_range(lp);
           init_state = rho_eg(:,lp,j1); 
           init_gg = V_eg_right_HEOM{j}*init_state; 
           init_ee = V_eg_left_HEOM{j}*init_state; 
           
            [t_out1,fun_out1] = ode45(rho_gg_prop,[0, tp_max],init_gg);
            [t_out2,fun_out2] = ode45(rho_gg_prop,[0,-tp_max],init_gg);
            t_out = [flip(t_out2),t_out1(2:end)];
            fun_out = [flip(fun_out1,2),fun_out2(:,2:end)];
            
            [t_out1,fun_out1] = ode45(rho_ee_prop,[0, tp_max],init_ee);
            [t_out2,fun_out2] = ode45(rho_ee_prop,[0,-tp_max],init_ee);
            t_out_ee = [flip(t_out2),t_out1(2:end)];
            fun_out_ee = [flip(fun_out1,2),fun_out2(:,2:end)];
            
            for om_lp1 = 1:length(omega_1)
            om1 = omega_1(om_lp1);
                for  om_lp2 = 1:length(omega_2)
                if ispumpprobe; om2 = om1; end %%PP
                
            to_int1 = repmat(conj(E2(t_out))*E1(t_out-t1)...
                        .*exp(1i*(om1-om_ex)*t1+1i*(om2-om1)*t_out) ...
                        ,size(fun_out,1)).* fun_out;
            int_cont = trapz(t_out,to_int1);
            
            to_int2 = repmat(conj(E2(t_out_ee))*E1(t_out_ee-t1)...
                        .*exp(1i*(om1-om_ex)*t1+1i*(om2-om1)*t_out_ee) ...
                        ,size(fun_out,1)).* fun_out_ee;
            int_cont2 = trapz(t_out_ee,to_int2);      
         if ispumpprobe         
rho_door_ee(:,om_lp1,j1,j2) =  rho_door_ee(:,om_lp1,j1,j2)+int_cont2;  
rho_door_gg(:,om_lp1,j1,j2) =  rho_door_gg(:,om_lp1,j1,j2)+int_cont;             
         else
rho_door_ee(:,om_lp1,lp2,j1,j2) = rho_door_ee(:,om_lp1,lp2,j1,j2)+int_cont2;  
rho_door_gg(:,om_lp1,lp2,j1,j2) = rho_door_gg(:,om_lp1,lp2,j1,j2)+int_cont;                       
         end
                end
            end
            
        end
        
       % else %adaptive 
       %     
       % end
    end
end

end
