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
 x=[1,0,0]; y=[0,1,0]; z=[0,0,1]; clear beam_param_set
if PP_setup
 % unique contributions from pump probe configuration
 beam_param_set{1} = cat(3,[x;x;x;x],[x;x;y;y]);
 F5=ori_F_calc(cat(3,[x;x;y;x;z],[y;y;y;x;z],[z;z;x;y;z])); 
  % Colinear:  xxxy (zzzz), xxyx (zzzz), Noncolinear:  zzyx (xxzz)
 beam_param_set{2} = cat(3,[x;x;y;x; z;z;z ],[ x;x;y;x ; y;y;z],... %
              [y;y;y;x; z;z;z],[ y;y;y;x; x;x;z],[z;z;x;y; x;x;z]); 
else
 % contributions from 2D configurations
  %beam_param_set{1} = cat(3,[x;x;y;y],[x;y;x;y],[x;y;y;x]);
  %take instead magic angle / cross peak / coherence
  ma = [1,sqrt(2),0]/sqrt(3); x1 = [cos(pi/3),sin(pi/3),0]; 
  x2 = [cos(pi/3),-sin(pi/3),0];
  beam_param_set{1} = cat(3,[x;x; ma; ma],[x1;x2;x;x],... %magic angle/ crosspeak
           [x;y;(x+y)/sqrt(2);(x-y)/sqrt(2)]); % coherence selec
  
 beam_param_set{2} = cat(3, [y;x;x;x; z;z;z], [x;y;x;x ;z;z;z], ...
  [x;x;y;x; z;z;z],[y;x;x;x; -z;-y;y], [x;y;x;x ;y;z;y], [x;x;y;x; y;y;z], ...
  [z;y;x;x; x;x;z],[x;z;y;x; -z;-x;x], [z;x;y;x; x;z;x]);   
  
  F5=ori_F_calc(cat(3, [y;x;x;x;z], [x;y;x;x;z], ...
  [x;x;y;x;z],[y;x;x;x;z], [x;y;x;x;z], [x;x;y;x;z], ...
  [z;y;x;x;x],[x;z;y;x;x], [z;x;y;x;x]));

%reduced set
%beam_param_set{1} = [x;x;x;x]; beam_param_set{2} = [y;x;x;x; z;z;z];
end
  F4=ori_F_calc(beam_param_set{1});

%%  Choose what to include etc
%Hamiltonian for sub sec from Ed's thesis
%DBVd [17113.9 319.374 30.4857 7.5811
%DBVc 319.374 17033.3 20.3238 -43.8736
%PCB158c 30.4857 20.3238 15807.4 3.3873
%MBVb 7.5811 -43.8736 3.3873 16372.0]
%From Eds
% 
% load('PC645_param_file.mat','H_site','mu','mdip','R',...
%     'mode_freq_and_HR_fact','mu_sd','pdm_shift')
% 
% %choose which sites I wish to include
% %MBVα18      PCBβ158     PCBβ82     DBVβ50/61  
% %MBVα18      PCBβ158     PCBβ82     DBVβ50/61   
% inc_sites = [false,false,false,true,false,false,false,true]; %which sites to include
% H_site = H_site(inc_sites,inc_sites); 
% mu = mu(inc_sites,:);  mdip = mdip(inc_sites,:);    R=R(inc_sites,:);
% mu_sd = mu_sd(inc_sites);  
% 
% N = sum(double(inc_sites));

%Just the dimer... I am not even sure these dipole moments and positions
%are the correct way round
%mu =[ -1.100 , -2.436 ,0.406 ;1.215 ,2.395,-0.056  ] ; mdip = mu*0; 
mu =[1.215 ,2.395,-0.056;  -1.100 , -2.436 ,0.406] ; mdip = mu*0; 
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
 mode_sel = [2,4]; %choose which modes to include in the spectral density
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
%t2_range_fs = linspace(0,1000,150);  
t2_range_fs = [0,eps];
 numvib = {[],[]};  
 Kappa = 1;  max_tier = 3;
 
if ~use_HEOM

save(params_for_wrapper,'Temperature_in_Kelvin','system_parameters_file','numvib',...
 'beam_param_set','pop_only','om1_rng','t2_range_fs','om3_rng','t1_coh','om2_coh',...
 'om3_coh','num_realisations')
else

save(params_for_wrapper,'Temperature_in_Kelvin','system_parameters_file','numvib',...
 'Kappa','max_tier','pop_only','beam_param_set','om1_rng','t2_range_fs','om3_rng',...
 't1_coh','om2_coh','om3_coh','num_realisations')    
    
end
%%
if use_HEOM
    twoD_with_HEOM_wrapper;
%save in a format that can be loaded easily
%save(file_name,'R1','R2','R3','R4','R5','R6','om1_rng','om3_rng','t2_range_fs','Drude_modes','BO_modes'...
%        ,'lin_spec','tmp_gg_t','tmp_ee_t','max_tier','Kappa','num_realisations','beam_param_set','chk','-v7.3');
else
    twoD_with_Redfield_wrapper;
%save(file_name,'R1','R2','R3','R4','R5','R6','om1_rng','om3_rng','t2_range_fs','BO_modes'...
%        ,'lin_spec','tmp_gg_t','tmp_ee_t','numvib','num_realisations','beam_param_set','chk','-v7.3');
  
end
%clear R1 R2 R3 R4 R5 R6 %this is too much data to keep in the RAM, load selectively

CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
CMRmap=interp1(1:9,CMRmap,1:0.25:9);
%xlabel('\tau (fs)','FontSize',16); ylabel('Probe frequency \omega_3 cm^{-1}','FontSize',14);
