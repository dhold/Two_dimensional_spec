%% Next bit

params_for_wrapper = 'Parameters/file_for_wrapper.mat';
Temperature_in_Kelvin = 297; 
system_parameters_file= 'Parameters/PC645_parameters.mat';
file_name = 'saved_data/2Dcd_PC645_2_site.mat'; %file name to save under

[t_scale, B, ~]= inv_cm_unit_sys(Temperature_in_Kelvin); 
t1_coh =[]; om2_coh=[]; om3_coh = []; %no coherence shit

%choose options about the system
pure_dephasing = false; %use pure dephasing rather than redfield or HEOM
use_HEOM = true; %if not selected will use redfield
inc_ud_modes = false; %include underdamped modes
PP_setup = false; %take pump probe configurations
pop_only = false; %only consider direct population transitions in exciton basis

%beam configurations which will be required

%%  Choose what to include etc

mu =[ -1.100 , -2.436 ,0.406 ;1.215 ,2.395,-0.056  ] ; mdip = mu*0; 

R = [20.417 , 2.139 ,  27.442;15.759 , -8.186 , 34.154];
H_site = [17113.9, 319.374; 319.374, 17033.3];
N=2; pdm_shift = []; mu_sd = mu*0; sd_mat = [0,0];

[PP, H_eig] = eig(H_site); EE = diag(H_eig);
lam_dru = [100,100]; gam_dru = [100,100];
 load('Parameters/PC645_param_file.mat', 'mode_freq_and_HR_fact')
 sDO = mode_freq_and_HR_fact(:,2); omega_DO = mode_freq_and_HR_fact(:,1);
  clear B0_modes  Drude_modes 
 for lp=1:N;  Drude_modes{lp} = [lam_dru(lp),gam_dru(lp)]; end
 
  gam_bro = 32; %damping
 mode_sel = [4,2]; %choose which modes to include in the spectral density
 nmodes = length(mode_sel); 
  %preassign empty
 for lp=1:N;  BO_modes{lp} = zeros(0,4); end
 for lp = mode_sel 
     BO_modes{1} = [BO_modes{1};[sDO(lp)*omega_DO(lp), gam_bro,omega_DO(lp),0]]; 
     BO_modes{2} = [BO_modes{2};[sDO(lp)*omega_DO(lp), gam_bro,omega_DO(lp),0]];
 end
     
 num_realisations = 1;
 if num_realisations == 1; sd_mat = sd_mat*0; end %no static disorder
 
save(system_parameters_file,'R','mu','H_site','BO_modes',...
        'Drude_modes','sd_mat','mdip','mu_sd','pdm_shift') 
%%
%range of frequencies to consider

om1_rng = {linspace(1.6e4,1.8e4,60),0}; %range in cm^(-1) of {freq,time}
om3_rng = linspace(1.6e4,1.8e4,200) ;
t2_range_fs = linspace(0,1000,150);  
%t2_range_fs = [0,eps];
 numvib = {[],[]};  
 Kappa = 1;  max_tier = 2;
 
bb = false; %don't calculate beyond dipole
save(params_for_wrapper,'Temperature_in_Kelvin','system_parameters_file','numvib',...
 'Kappa','max_tier','pop_only','om1_rng','t2_range_fs','om3_rng','bb',...
 't1_coh','om2_coh','om3_coh','num_realisations')    
%%

    twoD_w_HEOM_NOA_wrapper;

%clear R1 R2 R3 R4 R5 R6 %this is too much data to keep in the RAM, load selectively

%%
F4=ori_F_calc(cat(3,[x;x;x;x],[x;y;(x+y)/sqrt(2);(x-y)/sqrt(2)]));
 
t_pos = 1;
R_s = R2_s+R3_s+R5_s; %full rephasing spectra
R_nc = (R_s(:,t_pos,:,1,1)*F4(1,2) +  R_s(:,t_pos,:,2,1)*F4(2,2) +...
            R_s(:,t_pos,:,3,1)*F4(3,2))/30;

figure
contourf(om1_rng{1},om3_rng,squeeze(real(R_nc(:,1:end-1))))
%shading interp

%%

[~,tmp1] = min(abs(om1_rng{1}-1.75e4+0*200));
[~,tmp3] = min(abs(om3_rng-1.635e4));
k2 = 1;
R_nc = (R_s(tmp3,:,tmp1,1,1)*F4(1,k2) +  R_s(tmp3,:,tmp1,2,1)*F4(2,k2) +...
            R_s(tmp3,:,tmp1,3,1)*F4(3,k2))/30;
        
 figure; plot(t2_range,real(R_nc))       
 
 % Compute prony analysis
 N_order = 9; lg = t2_range_fs < 1000; %choose range
t_red= t2_range_fs(lg)*t_scale/1000; dt2 = (t_red(2)- t_red(1)); fs = 1/dt2;
a_list = zeros(N_order,1); spoles= a_list; 
tau_list = spoles; omega_list = spoles;
[t_scale, B, ~]= inv_cm_unit_sys(Temperature_in_Kelvin);


  imp_resp = real(R_nc(lg));

[num,den] = prony(imp_resp,N_order,N_order);
[r,p,k]= residuez(num,den); %find residues

a_list(:)= r(:); %amplitude list
spoles(:)= log(p(:))*fs; %poles
tau_list(:)= 1./imag(spoles(:));%decay factor
omega_list(:) = imag(spoles(:));
%signals are imag so will be made of cosines

figure
 curve = zeros(N_order,length(t_red));
    for lp = 1:N_order
%         if omega_list(lp) > eps
%             fct = 2;
%         elseif omega_list(lp) < -eps
%             fct = NaN;
%         else
%             fct = 1;
%         end
%        curve(lp,:) = fct*a_list(lp).*exp(spoles(lp)*t_red);   
curve(lp,:) = a_list(lp).*exp(spoles(lp)*t_red);   
    end
    plot(t_red,real(curve))
    xlabel('t_2 (fs)','FontSize',16); ylabel('ES signal, a.u.','FontSize',14);
    %for lp =1:N_order
    %    set('DisplayName',num2str(omega_list(lp)*1e3/t_scale))
    %end
 
%%
CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
CMRmap=interp1(1:9,CMRmap,1:0.25:9);
%xlabel('\tau (fs)','FontSize',16); ylabel('Probe frequency \omega_3 cm^{-1}','FontSize',14);
