function [rho_eg,rho_eg_ft]=rho_eg_t_HEOM2(N,t1_range,om_scale,...
                                    sup_op_ge,HL,H_trunc_1,tol)
%  Simplified version, works in exciton basis and assumes no explicit 
%   vibrations are included 
%  om_scale = typical frequency of exciton transfer
%  L_eg is -1i[H,rho_{eg}], plus extra acting on one specific tier of the
%  density matrix set of the HEOM.  But only the bit which mixes this block
%  onto itself obviously.  

if ~exist('tol','var')
    tol = [1e-8,1e-6]; %tolerances
end 
out_trunc = N*H_trunc_1; %truncate HEOM output at this tier
    

if iscell(t1_range) %pass cell to get the thing Fourier transformed
   
    om_1_rng = t1_range{2};  %range of carrier frequencies
    t1_coh =  t1_range{3};  %range of seperations to show at the end (can be empty)
    
    if ~isempty(om_1_rng)
    t1_range = t1_range{1}; dt = t1_range(2)-t1_range(1); 
    t_max = t1_range(end);  om1 = pi*(-1/dt:2/t_max:1/dt); %ft var
    %rho_eg_ft = zeros(N,1,H_trunc_1,length(om_1_rng),N);
    rho_eg_ft = zeros(N*H_trunc_1,length(om_1_rng),N);
    take_fft = true;
    else
    rho_eg_ft = []; take_fft = false;    
    end
    if ~isempty(t1_coh)
        %rho_eg = zeros(N*1,1,H_trunc_1,length(t1_coh),N);
        rho_eg = zeros(N*H_trunc_1,length(t1_coh),N);
    end
else
    %rho_eg = zeros(N,1,H_trunc_1,length(t1_range),N); 
    rho_eg = zeros(N*H_trunc_1,length(t1_range),N); 
    rho_eg_ft = [];  take_fft = false; t1_coh = [];
end


 
for j = 1:N %loop over sites

    om_ex = om_scale(j);%typical frequency of exciton transition
    
    sup_op_ge_sc = sup_op_ge + 1i*eye(N*HL)*om_ex; %scale with this freq
    rho_eg_prop = @(t,v) sup_op_ge_sc*v;   %equation to solve
    
rho_init = zeros(N*HL,1); rho_init(j) = 1;  %jth exciton transition
%solve the resulting ODE with the solver
rho_eg_tmp = OD_wrapper(t1_range,rho_eg_prop,rho_init,out_trunc,'ode45',tol);

if take_fft  %fft first variable for rephasing and nonrephasing signal
    rho_eg_tmp2 = fftshift(ifft(rho_eg_tmp,[],1),1);
    %interpolate to desired frequency range, don't extrapolate
    rng_to_interp = om_1_rng-om_ex;
    lg1 = min(om1)<rng_to_interp; lg2 = max(om1)>rng_to_interp; 
    rho_eg_tmp3 = zeros(N*HL,length(rng_to_interp));
    rho_eg_tmp3(:,lg1 & lg2) = interp1(om1,rho_eg_tmp2,rng_to_interp(lg1 & lg2),'pchip').';
    
    %reshape it into block matricies
   % rho_eg_ft(:,:,:,:,j) =  reshape(rho_eg_tmp2,[N,1,H_trunc_1,length(om_1_rng)]);
    rho_eg_ft(:,:,j) = rho_eg_tmp3;
end
if ~isempty(t1_coh)
    rho_eg_tmp  = interp1(t1_range,rho_eg_tmp ,t1_coh,'pchip');
   %rho_eg(:,:,:,:,j) =  reshape(rho_eg_tmp.',[N,1,H_trunc_1,length(t1_coh)]); 
   rho_eg(:,:,j) =  rho_eg_tmp.'; 
else %output full time range
 %rho_eg(:,:,:,:,j) =  reshape(rho_eg_tmp.',[N,1,H_trunc_1,length(t1_range)]);
    rho_eg(:,:,j) =  rho_eg_tmp.';
end
end