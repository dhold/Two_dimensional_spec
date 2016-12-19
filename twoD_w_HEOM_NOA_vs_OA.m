%this code is used to compare the results of twoD_with_HEOM and
%twoD_w_HEOM_NOA.  These two codes should be able to produce the same data
%for an orientation averaged system

params_for_wrapper = 'Parameters/file_for_wrapper.mat';
Temperature_in_Kelvin = 77; 
system_parameters_file= 'parameters/plenio_paper_parameters.mat';
file_name = 'saved_data/with_OA.mat'; %file name to save under
file_name2 = 'saved_data/without_OA.mat'; %file name to save under

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

  F4=ori_F_calc(beam_param_set{1});

%parameters for the system 
E_a = 12328;    E_b = 12472;   J = 70.7;    
H_site = [E_a,J;J,E_b];

[M_prj,E_ex]= eig(H_site);     E_ex = diag(E_ex); 
delta_Ex = E_ex(2)-E_ex(1); 

mu = [1, 0.5, 0;0, 1, 0];  
mu_ex = M_prj*mu;    mu_sd = [0,0]; 
R = 10^(-7)*[0, 1, 0;  0, 0, 1]; mdip = mu*0; pdm_shift = [];
sd_mat = 0.17*delta_Ex*[1,1];  %DOES contain static D

gamma_mode  = 0.005*delta_Ex;  om_mode = delta_Ex;
SHR = 0.02; om_0 = {om_mode,om_mode}; 

gamma_site = 0.025*delta_Ex;  %site based pure dephasing rate

%overdampeddamped modes
gam_dru = 100;
lam_dru = 10;

om1_rng = linspace(1.20e4,1.28e4,100); %range in cm^(-1) of {freq,time}
om3_rng = linspace(1.20e4,1.28e4,200) ;
t2_range_fs = linspace(0,20,3);  

 num_realisations = 1;
 if num_realisations == 1; sd_mat = sd_mat*0; end %no static disorder

 %convergence type parameters
 numvib = {3,3};
 Kappa = 2;  max_tier = 2;

    Drude_modes = {[lam_dru,gam_dru],[lam_dru,gam_dru]};


       % BO_modes{1} = [SHR*om_0{1}, gamma_mode,om_0{1},0]; 
       % BO_modes{2} = [SHR*om_0{2}, gamma_mode,om_0{2},0];
        BO_modes{1} = []; 
        BO_modes{2} = [];
% decided to calculate this stuff seperately 
%if PP_setup && 1==0;   freq_diff = 0;  beam_width = 0.05; %width in inv cm
%else;    freq_diff = []; beam_width = [];  end

save('parameters/plenio_paper_parameters.mat','R','mu','H_site','BO_modes',...
        'Drude_modes','sd_mat','mdip','mu_sd','pdm_shift')


save(params_for_wrapper,'Temperature_in_Kelvin','system_parameters_file','numvib',...
 'Kappa','max_tier','pop_only','beam_param_set','om1_rng','t2_range_fs','om3_rng',...
 't1_coh','om2_coh','om3_coh','num_realisations')    

%% WITH OA
    twoD_with_HEOM_wrapper;
%save in a format that can be loaded easily
save(file_name,'R1','R2','R3','R4','R5','R6','om1_rng','om3_rng','t2_range_fs','Drude_modes','BO_modes'...
        ,'lin_spec','tmp_gg_t','tmp_ee_t','max_tier','Kappa','num_realisations','beam_param_set','chk','-v7.3');

%% without OA

%params_for_wrapper = 'Parameters/file_for_wrapper2.mat';

%save(params_for_wrapper,'Temperature_in_Kelvin','system_parameters_file','numvib',...
% 'Kappa','max_tier','pop_only','beam_param_set','om1_rng','t2_range_fs','om3_rng',...
% 't1_coh','om2_coh','om3_coh','num_realisations') 

    
 twoD_w_HEOM_NOA_wrapper;
 save(file_name2,'R1_s','R2_s','R3_s','R4_s','R5_s','R6_s','om1_rng','om3_rng','t2_range_fs','Drude_modes','BO_modes'...
        ,'lin_spec','tmp_gg_t','tmp_ee_t','max_tier','Kappa','num_realisations','-v7.3');
    
%% Used to test if the exciton basis works ok, it does
k2 = 1; cmp_test = 4;
str = strcat('R',num2str(cmp_test),'_s'); 
R_s  = load('saved_data/without_OA_and_RR.mat',str); R_s = R_s.(str);
RR_s  = load('saved_data/without_OA3.mat',str); RR_s = RR_s.(str);
R_nc = (R_s(:,:,:,1,1)*F4(1,k2) +  R_s(:,:,:,2,1)*F4(2,k2) +  R_s(:,:,:,3,1)*F4(3,k2))/30;
RR_nc = (RR_s(:,:,:,1,1)*F4(1,k2) +  RR_s(:,:,:,2,1)*F4(2,k2) +  RR_s(:,:,:,3,1)*F4(3,k2))/30;

    toplot1 = squeeze(real(RR_nc(:,1,:)));
    toplot2 = squeeze(real(R_nc(:,1,:)));
    
%figure;     surf(om1_rng,om3_rng,toplot1)
%figure;     surf(om1_rng,om3_rng,toplot2)
figure;     surf(om1_rng,om3_rng,toplot2-toplot1)
   %%
% Test same for chial bits
R_c = R_nc*0; RR_c = RR_nc*0;
kb = 3;   tmp_beam = beam_param_set{2}; 
tmp_beam = tmp_beam(:,:,kb); tmp_beam = [tmp_beam;[0,0,1]]; 
for Lp = 1:4
    F_fct = tmp_beam(1:4,:); F_fct(Lp,:) = cross(F_fct(Lp,:),tmp_beam(4+Lp,:));
    F_fct = ori_F_calc(F_fct);
%R_c = R_c + (R_s(:,:,:,1,Lp+1)*F_fct(1) +  R_s(:,:,:,2,Lp+1)*F_fct(2) +  R_s(:,:,:,3,Lp+1)*F_fct(2))/30;
RR_c = RR_c +(RR_s(:,:,:,1,Lp+1)*F_fct(1) +  RR_s(:,:,:,2,Lp+1)*F_fct(2) +  RR_s(:,:,:,3,Lp+1)*F_fct(3))/30;
end
    toplot1 = squeeze(imag(RR_c(:,1,:)));
    
str = strcat('R',num2str(cmp_test)); RR = load(file_name,str); RR = RR.(str); 
RR = RR(:,:,:,3+k2);      toplot2 = squeeze(imag(RR(:,1,:)));
   % toplot2 = squeeze(imag(R_c(:,1,:)));
figure;     surf(om1_rng,om3_rng,toplot1)
figure;     surf(om1_rng,om3_rng,toplot2)

figure;     surf(om1_rng,om3_rng,toplot2/max(abs(toplot2(:)))-toplot1/max(abs(toplot1(:))))

%% Test polarizations    
    close all
k2 = 1; cmp_test = 4;
        
str = strcat('R',num2str(cmp_test)); RR = load(file_name,str); RR = RR.(str); RR = RR(:,:,:,k2); 
str = strcat('R',num2str(cmp_test),'_s'); R_s  = load(file_name2,str); R_s = R_s.(str);
R_nc = (R_s(:,:,:,1,1)*F4(1,k2) +  R_s(:,:,:,2,1)*F4(2,k2) +  R_s(:,:,:,3,1)*F4(3,k2))/30;
    toplot1 = squeeze(real(RR(:,1,:)));
    toplot2 = squeeze(real(R_nc(:,1,:)));
    
figure;     surf(om1_rng,om3_rng,toplot1/max(abs(toplot1(:))))
figure;     surf(om1_rng,om3_rng,toplot2/max(abs(toplot2(:))))
figure;     surf(om1_rng,om3_rng,toplot2/max(abs(toplot2(:)))-toplot1/max(abs(toplot1(:))))


%%
    R1_nc = (R1_s(:,:,:,1,1)*F4(1,k2) +  R1_s(:,:,:,2,1)*F4(2,k2) +  R1_s(:,:,:,3,1)*F4(3,k2))/30;
    R2_nc = (R2_s(:,:,:,1,1)*F4(1,k2) +  R2_s(:,:,:,2,1)*F4(2,k2) +  R2_s(:,:,:,3,1)*F4(3,k2))/30;
    R3_nc = (R3_s(:,:,:,1,1)*F4(1,k2) +  R3_s(:,:,:,2,1)*F4(2,k2) +  R3_s(:,:,:,3,1)*F4(3,k2))/30;
    R4_nc = (R4_s(:,:,:,1,1)*F4(1,k2) +  R4_s(:,:,:,2,1)*F4(2,k2) +  R4_s(:,:,:,3,1)*F4(3,k2))/30;
    R5_nc = (R5_s(:,:,:,1,1)*F4(1,k2) +  R5_s(:,:,:,2,1)*F4(2,k2) +  R5_s(:,:,:,3,1)*F4(3,k2))/30;
    R6_nc = (R6_s(:,:,:,1,1)*F4(1,k2) +  R6_s(:,:,:,2,1)*F4(2,k2) +  R6_s(:,:,:,3,1)*F4(3,k2))/30;
