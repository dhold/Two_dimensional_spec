function [V_ge,V_ef]=window_fun_HEOM(V_ex_eg,V_ex_fe,N,H_trunc,HL,tot_vib,t3_range,...
                       sup_op_ge,sup_op_ef, om_scale_ge,om_scale_ef,tol)
%  om_scale = typical frequency of exciton transfer
%  L_full contains all the processes acting on one specific tier of the
%  density matrix set of the HEOM.  But only the bit which mixes this block
%  onto itself obviously
% If HL =1 and H_mix_ge is just a matrix of zeros this is just a normal
% calculation with open quantum systems
%= -1i[H,rho] - sum_j Trunc_op_j [V_j,[V_j ,rho]] - n*gamma_n  rho_n


                         
if ~exist('tol','var')
    tol = [1e-7,1e-5]; %tolerances
end 

if iscell(t3_range) %pass cell to get the thing Fourier transformed
   
    om_3_rng = t3_range{2};  %range of carrier frequencies
    
    t3_range = t3_range{1}; dt = t3_range(2)-t3_range(1); 
    t_max = t3_range(end);  om3 = pi*(-1/dt:2/t_max:1/dt); %ft var
 
    %V_ge = zeros(N*tot_vib,tot_vib,H_trunc,length(om_3_range),N);
    V_ge = zeros(N*tot_vib^2*H_trunc,length(om_3_rng),N);
    %V_ef = zeros(N*(N-1)/2*tot_vib,tot_vib,H_trunc,length(om_3_range),N,N*(N-1)/2);
    V_ef = zeros(N^2*(N-1)/2*tot_vib^2*H_trunc,length(om_3_rng),N,N*(N-1)/2);

    take_fft = true;
else
  
   % V_ge = zeros(tot_vib,N*tot_vib,H_trunc,length(t3_range),N);
   V_ge = zeros(N*tot_vib^2*H_trunc,length(t3_range),N);  
   % V_ef = zeros(tot_vib,N*(N-1)/2*tot_vib,H_trunc,length(t3_range),N,N*(N-1)/2);
   V_ef = zeros(N^2*(N-1)/2*tot_vib^2*H_trunc,length(t3_range),N,N*(N-1)/2);

    take_fft = false; %calculate it in time space
end
 
for j = 1:N

    om_ex = om_scale_ge(j);%typical frequency of exciton transition
    
    sup_op_ge_sc = sup_op_ge - 1i*eye(N*tot_vib^2*HL)*om_ex; %scale with this freq
    V_ge_prop = @(t,v) mtimesx(sup_op_ge_sc,'n',v,'n');   
    
V_init = V_ex_eg(:,:,j); %jth exciton transition
%this is not trivial as the HEOM code works in the site basis 
out_trunc = numel(V_init)*H_trunc;  %save to truncation level
%reshape to column vector and pad with zeros to full size
V_init = reshape(V_init,numel(V_init),1);
V_init = [V_init;repmat(zeros(size(V_init)),[HL-1,1])]; 

%solve the resulting ODE with the solver
V_ge_tmp = (OD_wrapper(t3_range,V_ge_prop,V_init,out_trunc,'ode45',tol));
if take_fft
    V_ge_tmp2 = fftshift(ifft(V_ge_tmp,[],1),1);
    %interpolate to desired frequency range
    V_ge_tmp2 = interp1(om3+om_ex,V_ge_tmp2,om_3_rng,'pchip').';
    %reshape it into block matricies
    %V_ge(:,:,:,:,j) =  reshape(V_ge_tmp2,[tot_vib,N*tot_vib,H_trunc,length(om_3_range)]);
    V_ge(:,:,j) = V_ge_tmp2;
else
%reshape it into block matricies
   % V_ge(:,:,:,:,j) =  reshape(V_ge_tmp.',[tot_vib,N*tot_vib,H_trunc,length(t3_range)]);
    V_ge(:,:,j) =  V_ge_tmp.';
end
        for f = 1:N*(N-1)/2
            
    om_ex = om_scale_ef(j,f);%typical frequency of exciton transition
    
    sup_op_ef_sc = sup_op_ef - 1i*eye(N^2*(N-1)/2*tot_vib^2*HL)*om_ex; %scale with this freq
    V_ef_prop = @(t,v) mtimesx(sup_op_ef_sc,'n',v,'n');   
    
V_init = V_ex_fe(:,:,j,f); %jth exciton - fth double exciton transition
%this is not trivial as the HEOM code works in the site basis 
out_trunc = numel(V_init)*H_trunc;  %save to truncation level
%reshape to column vector and pad with zeros to full size
V_init = [reshape(V_init,numel(V_init),1);repmat(zeros(numel(V_init),1),[HL-1,1])];

%solve the resulting ODE with the solver
V_ef_tmp = (OD_wrapper(t3_range,V_ef_prop,V_init,out_trunc,'ode45',tol));

if take_fft
    
    V_ef_tmp2 = fftshift(ifft(V_ef_tmp,[],1),1);
    %interpolate to desired frequency range
    V_ef_tmp2 = interp1(om3+om_ex,V_ef_tmp2,om_3_rng,'pchip').';
    %reshape it into block matricies
    %V_ef(:,:,:,:,j,f) =  reshape(V_ef_tmp2,[tot_vib,N*(N-1)/2*tot_vib,H_trunc,length(om_3_range)]);
    V_ef(:,:,j,f) = V_ef_tmp2;
    
else
%reshape it into block matricies
   % V_ef(:,:,:,:,j,f) =  reshape(V_ef_tmp.',[tot_vib,N*(N-1)/2*tot_vib,H_trunc,length(t3_range)]);
    V_ef(:,:,j,f) =  V_ef_tmp.';
end
        end
end