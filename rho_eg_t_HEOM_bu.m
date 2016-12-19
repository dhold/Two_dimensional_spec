function [rho_eg,rho_eg_ft]=rho_eg_t_HEOM(V_ex_eg,rho_0,N,tot_vib,t1_range,om_scale,...
                                    sup_op_ge,HL,H_trunc_1,tol)
%  om_scale = typical frequency of exciton transfer
%  L_eg is -1i[H,rho_{eg}], plus extra acting on one specific tier of the
%  density matrix set of the HEOM.  But only the bit which mixes this block
%  onto itself obviously.  H_mix_op has the mixing a decay of higher tiers
%   - sum_j Trunc_op_j [V_j,[V_j ,rho]] - n*gamma_n  rho_n

if ~exist('tol','var')
    tol = [1e-7,1e-5]; %tolerances
end 

    
if iscell(t1_range) %pass cell to get the thing Fourier transformed
   
    om_1_rng = t1_range{2};  %range of carrier frequencies
    t1_coh =  t1_range{3};  %range of seperations to show at the end (can be empty)
    
    if ~isempty(om_1_rng)
    t1_range = t1_range{1}; dt = t1_range(2)-t1_range(1); 
    t_max = t1_range(end);  om1 = pi*(-1/dt:2/t_max:1/dt); %ft var
   % rho_eg_ft = zeros(N*tot_vib,tot_vib,H_trunc_1,length(om_1_rng),N);
   rho_eg_ft = zeros(N*tot_vib,tot_vib,H_trunc_1,length(om_1_rng),N);
    take_fft = true;
    else
    rho_eg_ft = []; take_fft = false;    
    end
    if ~isempty(t1_coh)
        rho_eg = zeros(N*tot_vib,tot_vib,H_trunc_1,length(t1_coh),N);
    end
else
   % rho_eg = zeros(N*tot_vib,tot_vib,H_trunc_1,length(t1_range),N); 
    rho_eg = zeros(N*tot_vib^2*H_trunc_1,length(t1_range),N); 
    rho_eg_ft = [];
    take_fft = false; t1_coh = [];
end
 
for j = 1:N

    om_ex = om_scale(j);%typical frequency of exciton transition
    
    sup_op_ge_sc = sup_op_ge+1i*eye(N*tot_vib^2*HL)*om_ex; %scale with this freq
    rho_eg_prop = @(t,v) sup_op_ge_sc*v;   %equation to solve
    
rho_init = V_ex_eg(:,:,j)*rho_0(1:tot_vib,1:tot_vib); %jth exciton transition
%this is not trivial as the HEOM code works in the site basis 
out_trunc = numel(rho_init)*H_trunc_1;  %save to truncation level
%reshape to column vector and pad with zeros to full size
rho_init = [reshape(rho_init,numel(rho_init),1);repmat(zeros(size(rho_init)),[HL-1,1])]; 

%solve the resulting ODE with the solver
rho_eg_tmp = (OD_wrapper(t1_range,rho_eg_prop,rho_init,out_trunc,'ode45',tol));

if take_fft  %fft first variable for rephasing and nonrephasing signal
    rho_eg_tmp2 = fftshift(ifft(rho_eg_tmp,[],1),1);
    %interpolate to desired frequency range
    rho_eg_tmp2 = interp1(om1+om_ex,rho_eg_tmp2,om_1_rng).';
    %reshape it into block matricies
    rho_eg_ft(:,:,:,:,j) =  reshape(rho_eg_tmp2,[N*tot_vib,tot_vib,H_trunc_1,length(om_1_rng)]);
end
if ~isempty(t1_coh)
    rho_eg_tmp  = interp1(t1_range,rho_eg_tmp ,t1_coh);
   rho_eg(:,:,:,:,j) =  reshape(rho_eg_tmp.',[N*tot_vib,tot_vib,H_trunc_1,length(t1_coh)]);  
else %output full time range
 rho_eg(:,:,:,:,j) =  reshape(rho_eg_tmp.',[N*tot_vib,tot_vib,H_trunc_1,length(t1_range)]);
end

%reshape it into block matricies
   

end