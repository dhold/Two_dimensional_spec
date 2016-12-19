function [V_gg,V_ee]=window_fun_HEOM3(N,H_trunc,HL,t3_range,...
                       sup_op_ge,sup_op_ef, om_scale_ge,om_scale_ef,tol)
                   warning('incomplete function')
 % simplified version that works in the exciton basis with no vibrations
%  om_scale = typical frequency of exciton transfer
%  L_full contains all the processes acting on one specific tier of the
%  density matrix set of the HEOM.  But only the bit which mixes this block
%  onto itself obviously
                
if ~exist('tol','var')
    tol = [1e-8,1e-6]; %tolerances
end 
out_trunc = N*H_trunc;  %save to truncation level
out_trunc2 = N^2*(N-1)/2*H_trunc;  %save to truncation level
if iscell(t3_range) %pass cell to get the thing Fourier transformed
   
    om_3_rng = t3_range{2};  %range of carrier frequencies
    
    t3_range = t3_range{1}; dt = t3_range(2)-t3_range(1); 
    t_max = t3_range(end);  om3 = pi*(-1/dt:2/t_max:1/dt); %ft var

    V_gg = zeros(H_trunc,length(om_3_rng),N,N);
    V_ee = zeros(N^2*H_trunc,length(om_3_rng),N,N*(N-1)/2,N,N*(N-1)/2);

    take_fft = true;
else

    V_gg = zeros(H_trunc,length(t3_range),N);
    V_ee = zeros(N^2*H_trunc,length(t3_range),N,1+N*(N-1)/2);

    take_fft = false; %calculate it in time space
end
 
for j = 1:N

    om_ex = om_scale_ge(j);%typical frequency of exciton transition
    
    sup_op_ge_sc = sup_op_ge + 1i*eye(N*HL)*om_ex; %scale with this freq
    V_ge_prop = @(t,v) mtimesx(sup_op_ge_sc,'n',v,'n');   
    
V_init = zeros(N*H_trunc,1); %jth exciton transition
V_init(j) = 1; %in exciton basis this is a simple initial condition

%solve the resulting ODE with the solver
V_ge_tmp = (OD_wrapper(t3_range,V_ge_prop,V_init,out_trunc,'ode45',tol));

    for j2 =1:N 

        %
        
        
    end
    
if take_fft
    V_ge_tmp2 = fftshift(ifft(V_ge_tmp,[],1),1);
    rng_to_interp = om_3_rng-om_ex;
    lg1 = min(om3)<rng_to_interp; lg2 = max(om3)>rng_to_interp; 
    V_ge_tmp3 = zeros(N*HL,length(rng_to_interp));
    %interpolate to desired frequency range, don't ffs extrapolate
    V_ge_tmp3(:,lg1&lg2) = interp1(om3,V_ge_tmp2,rng_to_interp(lg1&lg2),'pchip').';
    %reshape it into block matricies
    %V_ge(:,:,:,:,j) =  reshape(V_ge_tmp2,[1,N*1,H_trunc,length(om_3_range)]);
    V_ge(:,:,j) = V_ge_tmp3;
else
%reshape it into block matricies
   % V_ge(:,:,:,:,j) =  reshape(V_ge_tmp.',[1,N*1,H_trunc,length(t3_range)]);
    V_ge(:,:,j) =  V_ge_tmp.';
end
        for f = 1:N*(N-1)/2
            
    om_ex = om_scale_ef(j,f);%typical frequency of exciton transition
    
    sup_op_ef_sc = sup_op_ef + 1i*eye(N^2*(N-1)/2*HL)*om_ex; %scale with this freq
    V_ef_prop = @(t,v) mtimesx(sup_op_ef_sc,'n',v,'n');   
    
V_init = zeros(N,N*(N-1)/2); %jth exciton - fth double exciton transition
V_init(j,f) = 1;
%reshape to column vector and pad with zeros to full size
V_init = [reshape(V_init,[N^2*(N-1)/2,1]);zeros(N^2*(N-1)/2*(HL-1),1)];

%solve the resulting ODE with the solver
V_ef_tmp = (OD_wrapper(t3_range,V_ef_prop,V_init,out_trunc2,'ode45',tol));

if take_fft
    
    V_ef_tmp2 = fftshift(ifft(V_ef_tmp,[],1),1);
    rng_to_interp = om_3_rng-om_ex;
    lg1 = min(om3)<rng_to_interp; lg2 = max(om3)>rng_to_interp; 
    V_ef_tmp3 = zeros(N^2*(N-1)/2*HL,length(rng_to_interp));    
    %interpolate to desired frequency range
    V_ef_tmp3(:,lg1&lg2) = interp1(om3,V_ef_tmp2,rng_to_interp(lg1&lg2),'pchip').';

    %reshape it into block matricies
    %V_ef(:,:,:,:,j,f) =  reshape(V_ef_tmp2,[1,N*(N-1)/2*1,H_trunc,length(om_3_range)]);
    V_ef(:,:,j,f) = V_ef_tmp3;
    
else
%reshape it into block matricies
   % V_ef(:,:,:,:,j,f) =  reshape(V_ef_tmp.',[1,N*(N-1)/2*1,H_trunc,length(t3_range)]);
    V_ef(:,:,j,f) =  V_ef_tmp.';
end
        end
end