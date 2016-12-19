function [rho_fg]=rho_fg_HEOM (V_ex_fg,N,tot_vib,H_trunc,HL,rho_eg,...
                    sup_op_ge,t2_range,om2_rng,om_scale,simp_mix,tol ) 
     
  t1L = size(rho_eg,4); %length first signal in time / freq space 
  
  if isempty(om2_rng)
  t2L = length(t2_range);  %length to prop in t2     
  else
  t2L = length(om2_rng);  %length to prop in om2
    dt = t2_range(2) - t2_range(1); t_max= t2_range(end);
    om2 = pi*(-1/dt:2/t_max:1/dt);
  end
    out_trunc_fg = tot_vib^2*H_trunc*(N-1)*N/2;
    
    if simp_mix
       rho_fg =  zeros(tot_vib*(N-1)*N/2,tot_vib,H_trunc,t2L,t1L,N,(N-1)*N/2);
    else
       rho_fg =  zeros(tot_vib*(N-1)*N/2,tot_vib,H_trunc,t2L,t1L,N,N,(N-1)*N/2); 
    end
    
 for j1 = 1:N %loop over first interaction
     if simp_mix
         j2_rng = j1; %all others neglected, usually a safe approx!
     else
         j2_rng = 1:N;
     end
    for j2 = j2_rng %second interaction operator
        for f = 1:N*(N-1)/2

    om_sc = om_scale(j2,f);%typical frequency of exciton transition
    
    sup_op_ge_sc = sup_op_ge+1i*eye(N*tot_vib^2*HL)*om_sc; %scale with this freq
    rho_fg_prop = @(t,v) sup_op_ge_sc*v;   %equation to solve            
            
    %multiply by V_fg to get 0-2 coherence
    tmp_fg = mtimesx(V_ex_fg(:,:,j2,f),'C',rho_eg(:,:,:,:,j1));
    tmp_fg = reshape(tmp_fg,tot_vib^2*H_trunc_1,t1L);
    %loop over t1L, yes this WILL be painfully slow if this is long
    for lp = 1:length(t1L) 
        
        tmp_fg_sec = tmp_fg(:,lp);
        tmp_fg_t = (OD_wrapper(t2_range,rho_fg_prop,tmp_fg_sec,out_trunc_fg,'ode45',tol)).';
        if isempty(om2_rng)
        rho_fg(:,:,:,:,lp,j2,j1) = reshape(tmp_fg_t,[tot_vib*(N-1)*N/2,tot_vib,H_trunc_2,t2L]);
        else
        tmp_fg_ft = fftshift(ifft(tmp_fg_t,[],1),1);  
        tmp_fg_ft = interp1(om2+om_sc ,tmp_fg_ft,om2_rng);
        rho_fg(:,:,:,:,lp,j2,j1) = reshape(tmp_fg_ft,[tot_vib*(N-1)*N/2,tot_vib,H_trunc_2,t2L]);    
        end
    end
    
        end
    end
end       