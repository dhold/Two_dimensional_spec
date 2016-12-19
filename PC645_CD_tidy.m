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

  %small set
   beam_param_set{1} =[x;x;x;x]; beam_param_set{2} =[x;x;y;x; z;z;z];
%%  Choose what to include etc, incoherent rate from included to excluded
%{
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
%}

%This part deals the incoherent loss from the dimer subsystem to the rest
%of the system and ignores back coupling (i.e. uni-directional flow approx)
tic
HH = [17113.9 319.374 30.4857 7.5811
      319.374 17033.3 20.3238 -43.8736
      30.4857 20.3238 15807.4 3.3873
      7.5811 -43.8736 3.3873 16372.0]; %4 site Hamiltonian
lam_dru = [100,100]; gam_dru = [100,100];
 load('Parameters/PC645_param_file.mat', 'mode_freq_and_HR_fact')
  sDO = mode_freq_and_HR_fact(:,2); omega_DO = mode_freq_and_HR_fact(:,1);
  
  %Assume all 4 sites have same coloured bath
  gam_bro = 32; %decoherence rate of vib modes
  Drude_modes = [lam_dru,gam_dru]; 
  BO_modes = [sDO.*omega_DO, gam_bro.*(sDO.^0),omega_DO,sDO*0]; 
  %  BO_modes = [0,1,1,0]; %ignore these
  
sub_sys = [true,true,false,false]; N = sum(double(sub_sys));
[R_12,R_21,M_prj1] = Forster_outcouple(B,HH, sub_sys,Drude_modes,BO_modes,3000);
RR = sum(R_12,2); %rate at which pop leaves the exciton states
R_site = M_prj1'*diag(RR)*M_prj1; %relaxation tensor for sites

 % a Lindblad type operator but lossy
L_loss = -diag(RR)/2; %gives the same thing as before, just anti commutator
L_loss = kron(L_loss.',eye(length(RR)))+kron(eye(length(RR)),L_loss);
L_loss2 = -R_site/2;  %site bases decay terms, not actually diagonal
L_loss2 = kron(L_loss2.',eye(length(RR)))+kron(eye(length(RR)),L_loss2);

L_decay = -diag([0;RR;sum(RR)]); %include ground state and double excitons
L_decay = kron(L_decay.',eye(length(L_decay)))+kron(eye(length(L_decay)),L_decay);

% Full op that involves the extra state

%[sig_11,sig_21,sig_12,sig_22] = deal(HH*0); 
%sig_11(1,2+1) = 1; sig_12(1,2+2) = 1; sig_21(2,2+1) = 1; sig_22(2,2+2) = 1; 
%Lindblad = sig_11
%R_site2 = zeros(2,2);
%for k =1:2; for kk=1:2; for aa = 1:2; 
%R_site2(k,kk) =R_site2(k,kk) +  M_prj1(aa,k)*M_prj1(aa,kk)*RR(aa); end;end;end
toc
%% This bit is related to the dimer parameters
mu =[1.215 ,2.395,-0.056;  -1.100 , -2.436 ,0.406] ; mdip = mu*0; 
R = [20.417 , 2.139 ,  27.442;15.759 , -8.186 , 34.154];
H_site = [17113.9, 319.374; 319.374, 17033.3];
pdm_shift = []; mu_sd = mu*0; sd_mat = [0,0];

[PP, H_eig] = eig(H_site); EE = diag(H_eig);


  clear BO_modes  Drude_modes 
 for lp=1:N;  Drude_modes{lp} = [lam_dru(lp),gam_dru(lp)]; end
 
  gam_bro = 32; %damping / decoherence rate
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
        'Drude_modes','sd_mat','mdip','mu_sd','pdm_shift','L_decay') 
%%
%range of frequencies to consider

%om1_rng = linspace(1.6e4,1.8e4,60);%{linspace(1.6e4,1.8e4,60),0}; %range in cm^(-1) of {freq,time}
om1_rng = 8065.73*2.21 + linspace(-1e3,1e3,40);
om3_rng = 8065.73*2.01 + linspace(-1e3,1e3,90) ; %om3_rng = linspace(1.6e4,1.8e4,200) ;
Np2 = 100; t2_range_fs = linspace(0,1000,Np2);  
dt2 = (t2_range_fs(2)-t2_range_fs(1))/1000*t_scale; %t_in_ps *t_scale == t_in_inverse_cm
om2_rng = pi/dt2/(Np2-1)*((1-Np2):2:(Np2-1));
%t2_range_fs = [0,eps];
 numvib = {[],[]};  
 Kappa = 1;  max_tier = 2;
 
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
CMRmap_fine=interp1(1:9,CMRmap,1:0.25:9);
%xlabel('\tau (fs)','FontSize',16); ylabel('Probe frequency \omega_3 cm^{-1}','FontSize',14);

%% 
om2_rng2 = pi/dt2/(Np2-1)*(0:(Np2-1))*2;
tmp = R2(:,:,:,1)+R3(:,:,:,1)+R5(:,:,:,1);
tmp2 = fft(tmp,[],2); %tmp2 = fftshift(fft(tmp,[],2),2);
yy = sinc_filter(tmp2,dt2/pi*[700,950],[],2);
%y2 = ifft
filter_fun = exp(-(om2_rng-806.573).^2/2/350^2);
%filter_fun = ifftshift(filter_fun)/trapz(om2_rng,filter_fun);
tmp3 = tmp2.*repmat(filter_fun,[length(om3_rng),1,length(om1_rng)]);
tmp4 = ifft(tmp3,[],2); %tmp4 = ifftshift(ifft(tmp3,[],2),2);