
params_for_wrapper = 'Parameters/file_for_wrapper.mat';
Temperature_in_Kelvin = 77; 
[t_scale, B, ~]= inv_cm_unit_sys(Temperature_in_Kelvin); %time scaling and beta

t1_coh =[]; om2_coh=[]; om3_coh = []; %no coherence shit
clear beam_param_set

use_HEOM = true; 
inc_UD_mode  = false;

 % contributions to the nonchiral signal (not all indep)
 x=[1,0,0]; y=[0,1,0]; z=[0,0,1];
  F4=ori_F_calc(cat(3,[x;x;x;x],[x;x;y;y],[x;y;x;y],[x;y;y;x]));
 beam_param_set{1} = cat(3,[x;x;x;x],[x;x;y;y],[x;y;x;y],[x;y;y;x]);

  % contributions to the chiral signal (not all indep)
    F5=ori_F_calc(cat(3, [y;x;x;x;z], [x;y;x;x;z], ...
  [x;x;y;x;z],[y;x;x;x;z], [x;y;x;x;z], [x;x;y;x;z], ...
  [z;y;x;x;x],[x;z;y;x;x], [z;x;y;x;x],[x;x;x;-y;z],[x;x;x;-y;z]));
 beam_param_set{2} = cat(3, [y;x;x;x;z;z;z], [x;y;x;x;z;z;z], ...
  [x;x;y;x;z;z;z],[y;x;x;x;-z;-y;y], [x;y;x;x;y;z;y], [x;x;y;x;y;y;z], ...
  [z;y;x;x;x;x;z],[x;z;y;x;-z;-x;x], [z;x;y;x;x;z;x],...
  [z;z;z; x;x;x;-y],[y;y;z; x;x;x;-y]);                         


%same shit as before
E = 12000;  DE = 0; J = 100;
H_site = [E-DE/2,J;J,E+DE]; 
[M_prj,E_ex]= eig(H_site);   E_ex = diag(E_ex);  delta_Ex = E_ex(2)-E_ex(1);


mu = [1, 0, 0;0, 1, 0];   mu_ex = M_prj*mu;   
R = 10^(-7)*[0, 0, -1;0, 0, 1]; mdip = mu*0; pdm_shift = [];

gam_dru = 100;
lam_dru = 100; %important to change, can give diff bath each site
Drude_modes = {[lam_dru,gam_dru],[lam_dru,gam_dru]};

if inc_UD_mode

    mode_detune = 0;  
   om_0 = delta_Ex - mode_detune; 
  om_0 = {om_0_rng,om_0_rng}; 
  lam_bro = SHR*om_0_rng; 
 gam_bro =  pi*0.005*delta_Ex; 
 
 theta = asin(M_prj(2,1)); detune = om_0_rng-delta_Ex;
gtilde = sqrt(2)*cos(theta)*sin(theta)*sqrt(SHR*om_0_rng^2);
phi = atan(1/(-(detune)/2/gtilde+sqrt(1+(detune/2/gtilde)^2))); %mix angle   
Enew = (om_0_rng+delta_Ex)/2 + [-1,+1]*sqrt(gtilde^2+detune^2/4);
 
 numvib = {3,3};
%can take BO_modes out of Hamiltonian
if 1==0 %modes in Ham
 BO_modes{1} = [lam_bro, gam_bro,om_0_rng,numvib{1}];
 BO_modes{2} = [lam_bro, gam_bro,om_0_rng,numvib{2}];    
 fname = 'saved_data/2Dcd_homo_dimer_mode_in_Ham.mat' 
else %modes in spec density
BO_modes{1} = [lam_bro, gam_bro,om_0_rng,0]; 
BO_modes{2} = [lam_bro, gam_bro,om_0_rng,0];
 numvib = {0,0};
 fname = 'saved_data/2Dcd_homo_dimer_mode_in_SD.mat'
 end
else %no modes in system
BO_modes{1} = []; 
BO_modes{2} = [];
numvib = {0,0};
 fname = 'saved_data/2Dcd_homo_dimer_no_mode.mat' %#ok<*NASGU>
end

sd_mat = [0,0];%[0.17*delta_Ex,0.17*delta_Ex];  %Variance of site energies


%Range of frequencies to consider, second element is the time range
%Code calculates S(om_3, t_2, om_1 / t_1)
om1_rng = {E_ex,0};  %range of pump frequencies cm^(-1)
om3_rng = linspace(mean(E_ex)-3*delta_Ex,mean(E_ex)+3*delta_Ex,160) ;
 
t2_range_fs = linspace(0,2000,210);  
num_realisations = 1; %number of SD realisations to take, if no SD set to 1


save('parameters/plenio_paper_parameters.mat','R','mu','H_site','BO_modes',...
        'Drude_modes','sd_mat','mdip','mu_sd','pdm_shift')

Kappa =3;  max_tier = 4;  %use to achieve convergence

save(params_for_wrapper,'Temperature_in_Kelvin','system_parameters_file','numvib',...
 'Kappa','max_tier','pop_only','beam_param_set','om1_rng','t2_range_fs','om3_rng',...
 't1_coh','om2_coh','om3_coh','num_realisations')    
    

twoD_with_HEOM_wrapper; %call this function 
%save in a format that can be loaded easily
save(fname,'R1','R2','R3','R4','R5','R6','om1_rng','om3_rng'...
        ,'lin_spec','tmp_gg_t','tmp_ee_t','max_tier','Kappa',...,
        'num_realisations','beam_param_set','-v7.3');

clear R1 R2 R3 R4 R5 R6 %this is too much data to keep in the RAM, load selectively