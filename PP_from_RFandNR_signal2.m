
%ranges of the parameters
Temperature_in_Kelvin = 77;   
[t_scale, B, ~]= inv_cm_unit_sys(Temperature_in_Kelvin);
om1_rng = linspace(1.21e4,1.27e4,100); 
om3_rng = linspace(1.21e4,1.27e4,300);
t2_range_fs = linspace(0,2000,210);    t2_range = t2_range_fs*t_scale/1e3;
dt2 = t2_range_fs(2)-t2_range_fs(1); dt2_cm = dt2*t_scale/1e3;
fle = matfile('saved_data/plenio_w_exp_mode_indep_config_more_points2.mat');

tau_range_to_plot = t2_range(1:4:end);
om3_to_plot = om3_rng(1:10:end); om1_to_plot = om1_rng(1:10:end);

%Gaussian pulse parameters

sigma1_fs = 150; sigma1 = 150*t_scale/1e3; %pulse width in inverse cm
sigma2_fs = 100; sigma2 = 150*t_scale/1e3; 
max_sd = 10; %max Gaussian sds to take, 5 should probably be enough
num_pnts = 100; %number points in the t_range
if num_pnts ==0
t_range = -max_sd*sigma1:dt2_cm:max_sd*sigma1;
t_prime_range = -max_sd*sigma2:dt2_cm:max_sd*sigma2;
else
t_range = max_sd*sigma1*linspace(-1,1,num_pnts);
t_prime_range = max_sd*sigma2*linspace(-1,1,num_pnts);
end
[t_tilde,t_prime] = meshgrid(t_range,t_prime_range);

om_prime = pi*linspace(-1,1,length(t_prime_range))/(t_prime_range(2)-t_prime_range(1));
om = pi*linspace(-1,1,length(t_range))/(t_range(2)-t_range(1));
[om_p,omm] = meshgrid(om,om_prime);
%t_range = linspace(-7*max(sigma1,sigma2),7*max(sigma1,sigma2),300);
%write full possible range of t-t' values that we can have

E1 = exp(-t_prime.^2/2/sigma1^2)/sqrt(sqrt(pi)*sigma1);
E2 = exp(-t_tilde.^2/2/sigma2^2)/sqrt(sqrt(pi)*sigma2);
E1f = exp(-om_prime.^2*sigma1^2/2)/sqrt(sqrt(pi)/sigma1);
E2f = exp(-om.^2*sigma2^2/2)/sqrt(sqrt(pi)/sigma2);

%points to select
 om_pnts=[12320,12320;12320,12520;12520,12320;12520,12520];
om1_pnt = [1,1,2,2]; om3_pnt = [1,2,1,2]; 
scale_fct2 = max(max(max(abs(fle.R3(:,:,:,1)+fle.R4(:,:,:,1)))));
cmp = 3;

for lp = 1:3
%selectively load from the file
if lp==1
Sig = imag(fle.R3(1:length(om3_rng),:,1:length(om1_rng),cmp))/scale_fct2; 
Sig2 = imag(fle.R4(1:length(om3_rng),:,1:length(om1_rng),cmp))/scale_fct2; 
elseif lp==2
Sig = imag(fle.R1(1:length(om3_rng),:,1:length(om1_rng),cmp))/scale_fct2; 
Sig2 = imag(fle.R2(1:length(om3_rng),:,1:length(om1_rng),cmp))/scale_fct2; 
else
Sig = imag(fle.R5(1:length(om3_rng),:,1:length(om1_rng),cmp))/scale_fct2; 
Sig2 = imag(fle.R6(1:length(om3_rng),:,1:length(om1_rng),cmp))/scale_fct2; 
end
% for given values of "tau" need to take a 2D fft of 
% S(tau +t - t') E1(t) E2(t') 

%tmp = Sig(om3_sel,:,om1_sel); 

E_field_fact = repmat(E1.*E2,1,1,sz(2),sz(3));
%%
for om1_lp =1:length(om1_to_plot)
    om1 = om1(om1_to_plot);
    for om3_lp =1:length(om3_to_plot)
         om3 = om3(om3_to_plot);
        for lp1 = 1:length(om)
            om_part = om(lp1);
            phase_fct1 = exp(1i*t_range*om1);
            
            xlower1 = find(om1 +om_part > om1_rng,1,'last')-1;
            
            ylower1 = om3 +om_part - om3_rng(xlower1);
            yupper1 = om3_rng(xlower1+1)-(om3 + om_part);
            dy1 = om3_rng(xlower1+1)-om3_rng(xlower1);                
            
            for lp2 = 1:length(om_prime)
            omp = om_prime(lp1);    
            phase_fct2 = exp(1i*t_prime_range*om2);      
            
            xlower2 = find(om1 +om_part > om1_rng,1,'last')-1;
            
            ylower2 = om3 +om_part - om3_rng(xlower1);
            yupper2 = om3_rng(xlower1+1)-(om3 + om_part);
            dy2 = om3_rng(xlower1+1)-om3_rng(xlower1);             
           
            
for tau_lp = 1:length(tau_range_to_plot)
    
    %interpolate Sig to find the appropriate om1 +/- om' and om3+om
    Sig(
    
    tmp2 = zeros(length(t_range)*length(t_prime_range),sz(2),sz(3));
    tmpp2 = zeros(length(t_range)*length(t_prime_range),sz(2),sz(3));  
    
    tau_tmp = tau_range_to_plot(tau_lp); %value of tau in loop
    tt_rng = tau_tmp + t_tilde - t_prime; %range of t2   
    %Theta_mat = double(tt_rng >= 0); %Heaviside theta not actually needed
    tt_rng_flat = reshape(tt_rng,numel(tt_rng),1); %flatten
    [t_sorted,new_order ] = sort(tt_rng_flat); %sort
    t_lg = t_sorted >= 0 & t_sorted <= max(t2_range); %only actually use these
    
    tmp2(t_lg,:,:) = interpn(t2_range,om3_rng,om1_rng,tmp,...
                             t_sorted(t_lg),om3_rng,om1_rng);
    tmp2(~t_lg,:,:) = 0;        
    
    tmpp2(t_lg,:,:) = interpn(t2_range,om3_rng,om1_rng,tmpp,...
                             t_sorted(t_lg),om3_rng,om1_rng);
    tmpp2(~t_lg,:,:) = 0;         
    
    tmp2 = tmp2(new_order,:,:); % put these back in the correct order before sorting
    tmpp2 = tmpp2(new_order,:,:); 
    
    %now reshape to the correct size
    tmp2 = reshape(tmp2,length(t_tilde),length(t_prime),sz(2),sz(3));
    tmpp2 = reshape(tmpp2,length(t_tilde),length(t_prime),sz(2),sz(3));    
    %include the electric field contribution
    tmp2 = tmp2.*E_field_fact;  tmpp2 = tmpp2.*E_field_fact;
    phase_corre_fct = exp(-1i*om_p*t_prime_range(1)*(1-1/length(t_prime_range))...
                          -1i*omm*t_range(1)*(1-1/length(t_range)));
     phase_corre_fct =  repmat(phase_corre_fct,1,1,sz(2),sz(3));
    %phase_corre_fct = exp(1i*((length(xx)-1)/2:-1:(1-length(xx))/2 ...
    %                       + (length(yy)-1)/2:-1:(1-length(yy))/2));    
    tmp3 = fftshift(fft2(phase_corre_fct.*tmp2)); 
    tmpp3 = fftshift(fft2(phase_corre_fct.*tmpp2)); 
    %now I have to do the loop over the frequencies I actually require as
    %well sadly
    
    for om3_lp =1:length(om3_to_plot)
        om3 = om3_to_plot(om3_lp); 

      n=numel(om)-1; h=(om(end)-om(1))/length(om);
      tmp5 = zeros(size(tmp3,2),size(tmp3,4)); tmpp5=tmp5;
        for lp2 = 1:length(om) %perform integral using simpson
            om_part = om(lp2);
            xlower = find(om3 +om_part > om3_rng,1,'last')-1;
            
            ylower = om3 +om_part - om3_rng(xlower);
            yupper = om3_rng(xlower+1)-(om3 + om_part);
            dy = om3_rng(xlower+1)-om3_rng(xlower);
            tmp4 =  squeeze(E2f(lp2).*(tmp3(lp2,:,xlower,:)*ylower +...
                          tmp3(lp2,:,xlower+1,:)*yupper)/dy);
            tmpp4 =  squeeze(E2f(lp2).*(tmp3(lp2,:,xlower,:)*ylower +...
                          tmp3(lp2,:,xlower+1,:)*yupper)/dy);   
           if lp2 == 1 || lp2 == length(om)  
               tmp5 = tmp5 + squeeze(tmp4);
               tmpp5 = tmpp5 + squeeze(tmpp4);
           elseif 2*floor(lp2/2) ~= lp2
                tmp5 = tmp5 + 2*squeeze(tmp4);
                 tmpp5 = tmpp5 + 2*squeeze(tmpp4);
           else
               tmp5 = tmp5 + 4*squeeze(tmp4);
                tmpp5 = tmpp5 + 4*squeeze(tmpp4);
           end
        end
        tmp5 = tmp5*h/3;        tmpp5 = tmpp5*h/3;
        % now finally integrate over om_prime
        for om1_lp = 1:length(om1_to_plot)
        om1 = om1_to_plot(om1_lp); n=numel(om_prime)-1; 
        h=(om_prime(end)-om_prime(1))/length(om_prime);
        tmp6 = 0;
        for lp2 = 1:length(om_prime) %perform integral using simpson
            om_part = om_prime(lp2);
            xlower = find(om1 +om_part < om1_rng,'last');
            ylower = om1 +om_part - om1_rng(xlower);
            yupper = om1_rng(xlower+1)-(om1 + om_part);
            dy = om1_rng(xlower+1)-om1_rng(xlower);
            tmp4 =  squeeze(E1f(lp2).*(tmp5(lp2,xlower)*ylower +...
                          tmp5(lp2,xlower+1)*yupper)/dy);
            %the rephasing is at minus this value frequency        
            tmpp4 =  squeeze(E1f(length(om_prime)+1-lp2)...
                            .*(tmpp5(lp2,xlower)*ylower +...
                            tmpp5(lp2,xlower+1)*yupper)/dy);   
                        
           if lp2 == 1 || lp2 == length(om_prime)  
               tmp6 = tmp6 + squeeze(tmp4);
               tmpp6 = tmpp6 + squeeze(tmp4);
           elseif 2*floor(lp2/2) ~= lp2
                tmp6 = tmp6 + 2*squeeze(tmp4);
                tmpp6 = tmpp6 + 2*squeeze(tmp4);
           else
               tmp6 = tmp6 + 4*squeeze(tmp4);
               tmpp6 = tmpp6 + 4*squeeze(tmp4);
           end
        end
        tmp6 = tmp6*h/3; tmpp6 = tmpp6*h/3;
        PP_sig(om3_lp,tau_lp,om1_lp,lp) = tmp6+tmpp6; %FINALLY
        end        
    end
    toc
    

end
end