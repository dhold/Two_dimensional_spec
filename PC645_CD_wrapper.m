%% Next bit

params_for_wrapper = 'Parameters/file_for_wrapper.mat';
Temperature_in_Kelvin = 77; 
system_parameters_file= 'parameters/plenio_paper_parameters.mat';
file_name = 'saved_data/2Dcd_with_HEOM_PC645.mat'; %file name to save under

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
  ma = [1,sqrt(2),0]/sqrt(3); x1 = [cos(pi/3),sin(pi/3),0]; x2 = [cos(pi/3),-sin(pi/3),0];
  beam_param_set{1} = cat(3,[x;x; ma; ma],[x1;x2;x;x],... %magic angle/ crosspeak
           [x;y;(x+y)/sqrt(2);(x-y)/sqrt(2)]); % coherence selec
  
 beam_param_set{2} = cat(3, [y;x;x;x; z;z;z], [x;y;x;x ;z;z;z], ...
  [x;x;y;x; z;z;z],[y;x;x;x; -z;-y;y], [x;y;x;x ;y;z;y], [x;x;y;x; y;y;z], ...
  [z;y;x;x; x;x;z],[x;z;y;x; -z;-x;x], [z;x;y;x; x;z;x]);   
  
  F5=ori_F_calc(cat(3, [y;x;x;x;z], [x;y;x;x;z], ...
  [x;x;y;x;z],[y;x;x;x;z], [x;y;x;x;z], [x;x;y;x;z], ...
  [z;y;x;x;x],[x;z;y;x;x], [z;x;y;x;x]));
end
  F4=ori_F_calc(beam_param_set{1});

%parameters for the system 
E_a = 12328;    E_b = 12472;   J = 70.7;    
H_site = [E_a,J;J,E_b];

[M_prj,E_ex]= eig(H_site);     E_ex = diag(E_ex); 
delta_Ex = E_ex(2)-E_ex(1); 

mu = [1, 0.5, 0;0, 1, 0];   mu_ex = M_prj*mu;    mu_sd = [0,0]; 
R = 10^(-7)*[0, 1, 0;  0, 0, 1]; mdip = mu*0; pdm_shift = [];
sd_mat = 0.17*delta_Ex*[1,1];  %DOES contain static D

gamma_mode  = 0.005*delta_Ex;  om_mode = delta_Ex;
SHR = 0.02; om_0 = {om_mode,om_mode}; 

gamma_site = 0.025*delta_Ex;  %site based pure dephasing rate

%underdamped modes
gam_dru = 100;
lam_dru = gam_dru*B*gamma_site/2; %gives an unconvincingly low number :/
lam_dru = 30

%{
theta = asin(M_prj(2,1)); detune = om_mode-delta_Ex;
   gtilde = sqrt(2)*cos(theta)*sin(theta)*sqrt(SHR*om_mode^2);
   gprime = (cos(theta)^2-sin(theta)^2)*sqrt(SHR*om_mode^2);
 %lam_bro = coupg^2/om_mode = S_hr * om_mode
    phi = atan(1/(-(detune)/2/gtilde+sqrt(1+(detune/2/gtilde)^2))); %mix angle   

      Enew = (om_mode+delta_Ex)/2 + [-1,+1]*sqrt(gtilde^2+detune^2/4);

nmode = 3;
H_rel_vib = diag(om_mode*(1/2:(nmode+1/2))); M_prj2 = kron(M_prj,eye(nmode+1));
H_vib_ex = (diag(sqrt(1:nmode),1) + diag(sqrt(1:nmode),-1));    
H_vib_ex = gtilde*kron([1,0;0,-1],H_vib_ex) + gprime*kron([1,0;0,1],H_vib_ex);
H2 =  kron(diag(E_ex),eye(nmode+1)) +  kron(eye(2),H_rel_vib) + M_prj2'*H_vib_ex*M_prj2;
H2= (H2+H2')/2;  [a,b] = eig(H2); b=diag(b); 
    %}
%%
%range of frequencies to consider

om1_rng = {linspace(1.20e4,1.28e4,60),0}; %range in cm^(-1) of {freq,time}
om3_rng = linspace(1.20e4,1.28e4,200) ;
t2_range_fs = linspace(0,2000,210);  

 num_realisations = 1;
 if num_realisations == 1; sd_mat = sd_mat*0; end %no static disorder

 %convergence type parameters
 numvib = {3,3};
 Kappa = 2;  max_tier = 2;

if pure_dephasing 
    Drude_modes = {[gamma_site;0],[gamma_site;0]}; %pure dehphasing 
else
    Drude_modes = {[lam_dru,gam_dru],[lam_dru,gam_dru]};
end
%can take BO_modes out of Hamiltonian
%lam_bro = SHR*om_0_rng;  
if inc_ud_modes 
    if use_HEOM %no modes in Ham but in spec density
        %I could just include them in the Hamiltonian but I decided not to
        BO_modes{1} = [SHR*om_0{1}, gamma_mode,om_0{1},0]; 
        BO_modes{2} = [SHR*om_0{2}, gamma_mode,om_0{2},0];
    else  %modes in Ham
        BO_modes{1} = [SHR*om_0{1}, gamma_mode,om_0{1},numvib{1}];
        BO_modes{2} = [SHR*om_0{2}, gamma_mode,om_0{2},numvib{2}];   
    end
else
    BO_modes{1} = []; 
    BO_modes{2} = [];
end


% decided to calculate this stuff seperately 
%if PP_setup && 1==0;   freq_diff = 0;  beam_width = 0.05; %width in inv cm
%else;    freq_diff = []; beam_width = [];  end

save('parameters/plenio_paper_parameters.mat','R','mu','H_site','BO_modes',...
        'Drude_modes','sd_mat','mdip','mu_sd','pdm_shift')

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
clear R1 R2 R3 R4 R5 R6 %this is too much data to keep in the RAM, load selectively

CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
CMRmap=interp1(1:9,CMRmap,1:0.25:9);
%xlabel('\tau (fs)','FontSize',16); ylabel('Probe frequency \omega_3 cm^{-1}','FontSize',14);
