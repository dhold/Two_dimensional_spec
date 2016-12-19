%this version works out some parameters for CD signals as well
for jj =1:3
    
params_for_wrapper = 'Parameters/file_for_wrapper.mat';

system_parameters_file= 'parameters/plenio_paper_parameters.mat';
t1_coh =[]; om2_coh=[]; om3_coh = []; %no coherence shit

Temperature_in_Kelvin = 77;   
[t_scale, B, ~]= inv_cm_unit_sys(Temperature_in_Kelvin); %needed to convert parameters in 
%their equation into parameters for the HEOM

om1_rng = linspace(1.21e4,1.27e4,50); %range in cm^(-1)
t2_range = linspace(0,0.1586,100);
t2_range_fs = t2_range/t_scale*1e3;  %delay time t_2 in fs
om3_rng = linspace(1.21e4,1.27e4,300) ;


clear beam_param_set
pol1 = [1,0,0]; beam_param_set{1} = [pol1;pol1;pol1;pol1];
kvec1 = [0,1,1]/sqrt(2); kvec2 = [0,0,1];  pol2 = [1,1i,0]/sqrt(2);
beam_param_set{2} = [kvec1;kvec1;kvec2;pol1;pol1;pol2;conj(pol2)];
%beam_param_set =cat(3,[kvec1;kvec1;kvec2;pol1;pol1;pol1;conj(pol1)],...
%                      [kvec1;kvec1;kvec2;pol1;pol1;pol2;conj(pol2)]) ;

E_a = 12328;  E_b = 12472;  J = 70.7;
H_site = [E_a,J;J,E_b];

[M_prj,E_ex]= eig(H_site);
E_ex = diag(E_ex);  delta_Ex = E_ex(2)-E_ex(1);
%delta_Ex = delta_Ex*2*pi % is this where they went wrong / I did??

gamma_mode  = 0.005*delta_Ex;
SHR = 0.02; %Huang-Rhys factor for the mode
om_0_rng = [0.75,1,1.25]*delta_Ex; %om_0 which was considered

%calculate splitting of first two excited modes from H_JC
detune = om_0_rng(jj)-delta_Ex; 
 gtilde = sqrt(2)*prod(M_prj(:,2))*sqrt(SHR)*om_0_rng(jj);
    phi = atan(1/(-detune/2/gtilde+sqrt(1+(detune/2/gtilde)^2)));   
   Enew = (om_0_rng(jj)+delta_Ex)/2 + [-1,+1]*sqrt(gtilde^2+detune^2/4);
      
  om_0 = {om_0_rng(jj),om_0_rng(jj)};

 lam_bro = SHR*om_0_rng(jj);
 lam_bro = lam_bro*pi
 
 gam_bro = gamma_mode;
 BO_modes{1} = [lam_bro, gam_bro,om_0_rng(jj),0];
 BO_modes{2} = [lam_bro, gam_bro,om_0_rng(jj),0];

gamma_site = 0.025*delta_Ex; %pure dephasing for each site
%relate the pure dephasing for each site to a Drude spectral density
%2 \lambda / \beta \gamma, assume \gamma = 100 and thus
%  \lambda = \beta \gamma gamma_site / 2
if 1 ==1
gam_dru = 100;
lam_dru = gam_dru*B*gamma_site/2; %gives an unconvincingly low number :/
lam_dru = lam_dru*2*pi %increase by factor of 2*pi
Drude_modes{1} = [lam_dru,gam_dru];
Drude_modes{2} = [lam_dru,gam_dru];
else %just include some pure dephasing shit
    %passing in this format will be pure dephasing
    Drude_modes{1} = [gamma_site*2*pi;0]
    Drude_modes{2} = [gamma_site*2*pi;0]
end

mu = [1, 0.5, 0;0, 1, 0];
mu_ex = M_prj*mu;
mu_sd = [0,0]; sd_mat = [0,0]; %no static disorder parameters
num_realisations = 1; numvib = {3,3};
R = 10^(-7)*[0, 0, 1;   0, 0, 0]; mdip = mu*0; %beyond dipole stuff

pdm_shift = []; %no shift to double excited state
%% Save
save(params_for_wrapper,'Temperature_in_Kelvin','system_parameters_file',...
 'beam_param_set','om1_rng','t2_range_fs','om3_rng','t1_coh','om2_coh','om3_coh','numvib','num_realisations') 

save('parameters/plenio_paper_parameters.mat','R','mu','H_site','BO_modes',...
        'Drude_modes','sd_mat','mdip','mu_sd','pdm_shift')

    %% run the code
    t_verbose = true;
    %twoD_with_HEOM_wrapper;
    twoD_with_Redfield_wrapper;
    %% save the response functions
    
 %   save(strcat('saved_data/plenio_w_heom_paramset_',num2str(jj),'.mat'),...
  %      'R1','R2','R3','R4','R5','R6','lin_spec','tmp_gg_t','tmp_ee_t','max_tier','Kappa',...
  %      'Drude_modes','BO_modes')
    save(strcat('saved_data/plenio5_wcd_Redf_paramset_',num2str(jj),'.mat'),...
        'R1','R2','R3','R4','R5','R6','lin_spec','tmp_gg_t','tmp_ee_t','numvib',...
    'Drude_modes','BO_modes','opts');    
    
end
%%
%figure
 for jj =1:3 %plot the graphs from plenios paper
    
    load(strcat('saved_data/plenio5_wcd_Redf_paramset_',num2str(jj),'.mat'),'R3')
cmp=2;
scale_fct2 = max(max(max(abs(real(R3(:,:,:,cmp))))));
pnt = 1%length(t2_range);
%subplot(3,1,jj)
figure
pcolor(om3_rng,om1_rng,(real(squeeze(R3(:,pnt,:,cmp)))).'/scale_fct2)
shading interp
xlabel('\omega_3 (cm^{-1})')
ylabel('\omega_1 (cm^{-1})')
 end
%%
 [a,pnt2] = min(abs(t2_range - 1.6/delta_Ex)); %find min 
load(strcat('saved_data/plenio5_wcd_Redf_paramset_',num2str(2),'.mat'),'R3','R2','R5')
 %%
pnt2=pnt2+10 %iterate by 1 to see tdep, pretty cool!
cmp2 = 2;
%load(strcat('test_data.mat'),'R3','R2','R5')
figure
pcolor(om3_rng,om1_rng,squeeze((real(R3(:,pnt2,:,cmp2)+R2(:,pnt2,:,cmp2)+R5(:,pnt2,:,cmp2)))).'/scale_fct2)
shading interp
xlabel('\omega_3')
ylabel('\omega_1')

%% also plot the slices
om_pnts=[12320,12320;12320,12520;12520,12320;12520,12520]; %om1,om3
cmp2  =2;
%find closest points to these
to_plot1 = zeros(size(om_pnts,1),length(t2_range)); 
to_plot2 = to_plot1; to_plot3 = to_plot2;
for lp=1:4
    [~,om1_pnt]=min(abs(om_pnts(lp,1)-om1_rng));
    [~,om3_pnt]=min(abs(om_pnts(lp,2)-om3_rng));
    to_plot1(lp,:) = real(R3(om3_pnt,:,om1_pnt,cmp2))/scale_fct2;
    to_plot2(lp,:) = real(R2(om3_pnt,:,om1_pnt,cmp2))/scale_fct2;
    to_plot3(lp,:) = real(R5(om3_pnt,:,om1_pnt,cmp2))/scale_fct2;                 
                 
end
%%
figure
plot(t2_range/t_scale*1000,to_plot1)
xlabel('t_2 (fs)')

%%
figure
plot(t2_range/t_scale*1000,to_plot2+to_plot3)
xlabel('t_2 (fs)')
%%
figure
plot(t2_range/t_scale*1000,to_plot2)
xlabel('t_2 (fs)')
figure
plot(t2_range/t_scale*1000,to_plot3)
xlabel('t_2 (fs)')
%%
if 1==1
%% Next bit

params_for_wrapper = 'Parameters/file_for_wrapper.mat';
Temperature_in_Kelvin = 77; system_parameters_file= 'parameters/plenio_paper_parameters.mat';
[t_scale, B, ~]= inv_cm_unit_sys(Temperature_in_Kelvin); 
t1_coh =[]; om2_coh=[]; om3_coh = []; %no coherence shit
clear beam_param_set
%{
pol1 = [1,0,0]; beam_param_set{1} = [pol1;pol1;pol1;pol1];
kvec1 = [0,1,1]/sqrt(2); kvec2 = [0,0,1];  pol2 = [1,1i,0]/sqrt(2); 
theta_beam = atan(sqrt(2)); kvec11 = [0,cos(theta_beam),sin(theta_beam)];
beam_param_set{2} = cat(3,[kvec1;kvec1;kvec2;pol1;pol1;pol2;conj(pol2)],... %
                          [kvec1;kvec1;kvec2;pol1;pol1;conj(pol2);pol2],... % pointless
                          [kvec1;kvec1;kvec2;pol1;pol1;pol1;[0,1,0]],...
                          [kvec1;kvec1;kvec2;pol1;pol1;[0,1,0];pol1],...
                          [kvec11;kvec11;kvec2;pol1;pol1;pol1;[0,1,0]],...
                         [kvec11;kvec11;kvec2;pol1;pol1;[0,1,0];pol1]);
%}
%use_HEOM = true; 
use_HEOM = false; 
PP_setup = false;
 x=[1,0,0]; y=[0,1,0]; z=[0,0,1];
if PP_setup
 % unique contributions from pump probe configuration

  F4=ori_F_calc(cat(3,[x;x;x;x],[x;x;y;y]));
 beam_param_set{1} = cat(3,[x;x;x;x],[x;x;y;y]);
 F5=ori_F_calc(cat(3,[x;x;y;x;z],[y;y;y;x;z],[z;z;x;y;z]));       
 beam_param_set{2} = cat(3,[x;x;y;x; z;z;z ],[ x;x;y;x ; y;y;z],... %
              [y;y;y;x; z;z;z],[ y;y;y;x; x;x;z],[z;z;x;y; x;x;z]); 
else
 % unique contributions from 2D configurations
  %beam_param_set{1} = cat(3,[x;x;y;y],[x;y;x;y],[x;y;y;x]);
  %take instead magic angle / cross peak / coherence
  ma = [1,sqrt(2),0]/sqrt(3); x1 = [cos(pi/3),sin(pi/3),0]; x2 = [cos(pi/3),-sin(pi/3),0];
  beam_param_set{1} = cat(3,[x;x; ma; ma],[x1;x2;x;x],... %magic angle/ crosspeak
           [x;y;(x+y)/sqrt(2);(x-y)/sqrt(2)]); % coherence selec
  F4=ori_F_calc(beam_param_set{1});
  
 beam_param_set{2} = cat(3, [y;x;x;x; z;z;z], [x;y;x;x ;z;z;z], ...
  [x;x;y;x; z;z;z],[y;x;x;x; -z;-y;y], [x;y;x;x ;y;z;y], [x;x;y;x; y;y;z], ...
  [z;y;x;x; x;x;z],[x;z;y;x; -z;-x;x], [z;x;y;x; x;z;x]);   
  
  F5=ori_F_calc(cat(3, [y;x;x;x;z], [x;y;x;x;z], ...
  [x;x;y;x;z],[y;x;x;x;z], [x;y;x;x;z], [x;x;y;x;z], ...
  [z;y;x;x;x],[x;z;y;x;x], [z;x;y;x;x]));
end                    
 % Colinear:  xxxy (zzzz), xxyx (zzzz),
 % Noncolinear:  zzyx (xxzz)
 %  (adding y to dir of pump (changing angle just drops signal amplitude)
 
%same shit as before
E_a = 12328;E_b = 12472;J = 70.7;H_site = [E_a,J;J,E_b];
[M_prj,E_ex]= eig(H_site);E_ex = diag(E_ex);  delta_Ex = E_ex(2)-E_ex(1);
gamma_mode  = pi*0.005*delta_Ex; SHR = 0.02; om_0_rng = delta_Ex; 
  om_0 = {om_0_rng,om_0_rng}; 
  lam_bro = SHR*om_0_rng;   gam_bro = gamma_mode; 

gamma_site = 0.025*delta_Ex; 
Drude_modes = {[gamma_site*2*pi;0],[gamma_site*2*pi;0]}; %pure dehphasing
mu = [1, 0.5, 0;0, 1, 0];mu_ex = M_prj*mu;    mu_sd = [0,0]; 
R = 10^(-7)*[0, 1, 0;0, 0, 1]; mdip = mu*0; pdm_shift = [];

sd_mat = [0.17*delta_Ex,0.17*delta_Ex];  %DOES contain static D
sd_mat = sd_mat*0 % doesn't if you uncomment this
%{
theta = asin(M_prj(2,1)); detune = om_0_rng-delta_Ex;
   gtilde = sqrt(2)*cos(theta)*sin(theta)*sqrt(SHR*om_0_rng^2);
   gprime = (cos(theta)^2-sin(theta)^2)*sqrt(SHR*om_0_rng^2);
 %lam_bro = coupg^2/omegavib = S_hr * omegavib
    phi = atan(1/(-(detune)/2/gtilde+sqrt(1+(detune/2/gtilde)^2))); %mix angle   

      Enew = (om_0_rng+delta_Ex)/2 + [-1,+1]*sqrt(gtilde^2+detune^2/4);

%%
H_rel_vib = diag(om_0_rng*(1/2:(7+1/2))); M_prj2 = kron(M_prj,eye(8));
H_vib_ex = (diag(sqrt(1:7),1) + diag(sqrt(1:7),-1));    
H_vib_ex = gtilde*kron([1,0;0,-1],H_vib_ex) + gprime*kron([1,0;0,1],H_vib_ex);
H2 =  kron(diag(E_ex),eye(8)) +  kron(eye(2),H_rel_vib) + M_prj2'*H_vib_ex*M_prj2;
H2= (H2+H2')/2;  [a,b] = eig(H2); b=diag(b); 
      %}
%%
%much smaller range of frequencies to consider
%om1_rng = [12320,12520]; om3_rng = om1_rng; 
%om1_range = [12320,12520]; t1_rng = 0;  om1_rng = {om1_range,t1_rng};

%om3_rng = linspace(1.21e4,1.27e4,300);

%Choose points distributed as a Lorentzian about the peaks?
om1_rng = linspace(1.20e4,1.28e4,80); %range in cm^(-1)
om3_rng = linspace(1.20e4,1.28e4,200) ;
% om1_rng = {[],0};

if 1==1 %use Redfield or HEOM
gam_dru = 100;
lam_dru = gam_dru*B*gamma_site/2; %gives an unconvincingly low number :/
lam_dru = lam_dru*2*pi %increase by factor of 2*pi
Drude_modes{1} = [lam_dru,gam_dru];
Drude_modes{2} = [lam_dru,gam_dru];
end

t2_range_fs = linspace(0,2000,210);   num_realisations = 1;
 numvib = {3,3};
%can take BO_modes out of Hamiltonian
%lam_bro=lam_bro*pi
%gam_bro=gam_bro*pi
if ~use_HEOM %modes in Ham
 BO_modes{1} = [lam_bro, gam_bro,om_0_rng,numvib{1}];
 BO_modes{2} = [lam_bro, gam_bro,om_0_rng,numvib{2}];    
elseif 1==1 %no modes in Ham but in spec density
BO_modes{1} = [lam_bro, gam_bro,om_0_rng,0]; 
BO_modes{2} = [lam_bro, gam_bro,om_0_rng,0];
else %no modes
BO_modes{1} = []; 
BO_modes{2} = [];
end


if PP_setup && 1==0;   freq_diff = 0;  beam_width = 0.05; %width in inv cm
else;    freq_diff = []; beam_width = [];  end

save('parameters/plenio_paper_parameters.mat','R','mu','H_site','BO_modes',...
        'Drude_modes','sd_mat','mdip','mu_sd','pdm_shift')
    pop_only = false
    
if ~use_HEOM


save(params_for_wrapper,'Temperature_in_Kelvin','system_parameters_file','numvib',...
 'beam_param_set','pop_only','om1_rng','t2_range_fs','om3_rng','t1_coh','om2_coh',...
 'om3_coh','num_realisations','freq_diff','beam_width')
else
    Kappa =3;  max_tier = 4; 

save(params_for_wrapper,'Temperature_in_Kelvin','system_parameters_file','numvib',...
 'Kappa','max_tier','pop_only','beam_param_set','om1_rng','t2_range_fs','om3_rng',...
 't1_coh','om2_coh','om3_coh','num_realisations','freq_diff','beam_width')    
    
end
%%
if use_HEOM
twoD_with_HEOM_wrapper;
%save in a format that can be loaded easily
save('saved_data/PPcd_with_HEOM_with_mode.mat','R1','R2','R3','R4','R5','R6','om1_rng','om3_rng'...
        ,'lin_spec','tmp_gg_t','tmp_ee_t','max_tier','Kappa','num_realisations','beam_param_set','t2_range_fs','chk','-v7.3');
else
    twoD_with_Redfield_wrapper;
  +
  
enj1?d
clear R1 R2 R3 R4 R5 R6 %this is too much data to keep in the RAM, load selectively
%%
%  load('saved_data/plenio_w_exp_mode_indep_config_more_points2.mat','R1','R2','R3','R4','R5','R6'...
%          ,'lin_spec','tmp_gg_t','tmp_ee_t');
      
%fle = matfile('saved_data/plenio_w_exp_mode_indep_config_more_points2.mat');
% too much data to load all at once really so just load what you need to
% plot
fle = matfile('saved_data/PPcd_with_HEOM_with_mode.mat');
%fle = matfile('saved_data/plenio_cd_config_with_RF_nomode.mat')
 %%     
 om_pnts=[12320,12320;12320,12520;12520,12320;12520,12520];
om1_pnt = [1,1,2,2]; om3_pnt = [1,2,1,2]; 
scale_fct2 = max(max(max(abs(fle.R3(:,:,:,1)+fle.R4(:,:,:,1)))));
for cmp = (size(beam_param_set{1},3)+(1:size(beam_param_set{2},3)))
    %%
%scale_fct2 = max(max(max(abs(R4(:,:,:,1)+R3(:,:,:,1)))));
to_plot1 = zeros(size(om_pnts,1),length(t2_range_fs)); 
to_plot2 = to_plot1; to_plot3 = to_plot2;
tau = inf;%1e-2; %sd of pump pulse
GSB = imag(fle.R3(:,:,:,cmp)+fle.R4(:,:,:,cmp))/scale_fct2; 
SE = imag(fle.R1(:,:,:,cmp)+fle.R2(:,:,:,cmp))/scale_fct2; 
ESA = imag(fle.R5(:,:,:,cmp)+fle.R6(:,:,:,cmp))/scale_fct2;
%%
if iscell(om1_rng)
    om1_range = om1_rng{1};
    t1_range = om1_rng{2};
else
    om1_range = om1_rng;
end

for lp=1:4
    [~,om1_pnt]=min(abs(om1_range-om_pnts(lp,1)));
    [~,om3_pnt]=min(abs(om3_rng-om_pnts(lp,2)));
    
GSB_cmp = GSB(om3_pnt,:,om1_pnt);
SE_cmp = SE(om3_pnt,:,om1_pnt);
ESA_cmp = ESA(om3_pnt,:,om1_pnt);
    if tau == inf %infinite -> frequency resolved
    to_plot1(lp,:) = GSB_cmp;
    to_plot2(lp,:) = SE_cmp;
    to_plot3(lp,:) = ESA_cmp; 
    elseif tau == 0 %0 -> time resolved |E(omega)|^2 ~ const
        if iscell(om1_rng)
            lg=t1_range == 0 ;
    to_plot1(lp,:) = GSB(  om3_pnt,:,length(om1_rng)+1);  
    to_plot2(lp,:) = SE(  om3_pnt,:,length(om1_rng)+1);  
    to_plot3(lp,:) = ESA(  om3_pnt,:,length(om1_rng)+1);          
        else
    to_plot1(lp,:) = trapz(om1_rng,GSB_cmp,3);     
    to_plot2(lp,:) = trapz(om1_rng,SE_cmp,3);    
    to_plot3(lp,:) = trapz(om1_rng,ESA_cmp,3);
        end
    else %could also do this with time range if I wanted
        om_pump = om1_rng(om1_pnt);
    to_plot1(lp,:) = PP_from_PE(permute(GSB_cmp,[3,2,1,4])...
                    ,tau,om_pump,om1_rng,[],[],[]);
    to_plot2(lp,:) = PP_from_PE(permute(SE_cmp,[3,2,1,4])...
                    ,tau,om_pump,om1_rng,[],[],[]);
    to_plot3(lp,:) = PP_from_PE(permute(ESA_cmp,[3,2,1,4])...
                    ,tau,om_pump,om1_rng,[],[],[]);   
    end      
    %{
    if tau == inf %infinite -> frequency resolved
    to_plot1(lp,:) = real(R3(om3_pnt,:,om1_pnt,cmp)+...
                          R4(om3_pnt,:,om1_pnt,cmp))/scale_fct2;
    to_plot2(lp,:) = real(R2(om3_pnt,:,om1_pnt,cmp)+...
                          R1(om3_pnt,:,om1_pnt,cmp))/scale_fct2;
    to_plot3(lp,:) = real(R5(om3_pnt,:,om1_pnt,cmp)+...
                          R6(om3_pnt,:,om1_pnt,cmp))/scale_fct2; 
    elseif tau == 0 %0 -> time resolved |E(omega)|^2 ~ const
    to_plot1(lp,:) = trapz(om1_rng,real(R3(om3_pnt,:,:,cmp)+...
                          R4(om3_pnt,:,:,cmp))/scale_fct2,3);     
    to_plot2(lp,:) = trapz(om1_rng,real(R2(om3_pnt,:,:,cmp)+...
                          R1(om3_pnt,:,:,cmp))/scale_fct2,3);    
    to_plot3(lp,:) = trapz(om1_rng,real(R5(om3_pnt,:,:,cmp)+...
                          R6(om3_pnt,:,:,cmp))/scale_fct2,3);                         
    else
        om_pump = om1_rng(om1_pnt);
    to_plot1(lp,:) = PP_from_PE(permute(real(R3(om3_pnt,:,:,cmp)+ R4(om3_pnt,:,:,cmp))...
                    /scale_fct2,[3,2,1,4]),tau,om_pump,om1_rng,[],[],[]);
    to_plot2(lp,:) = PP_from_PE(permute(real(R1(om3_pnt,:,:,cmp)+ R2(om3_pnt,:,:,cmp))...
                    /scale_fct2,[3,2,1,4]),tau,om_pump,om1_rng,[],[],[]);
    to_plot3(lp,:) = PP_from_PE(permute(real(R5(om3_pnt,:,:,cmp)+ R6(om3_pnt,:,:,cmp))...
                    /scale_fct2,[3,2,1,4]),tau,om_pump,om1_rng,[],[],[]);   
    end  
    %}
end
to_plot_ground{cmp} = to_plot1;
to_plot_excited{cmp} = to_plot2+to_plot3;
end
CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
CMRmap=interp1(1:9,CMRmap,1:0.25:9);
%xlabel('\tau (fs)','FontSize',16); ylabel('Probe frequency \omega_3 cm^{-1}','FontSize',14);

%%
figure1 = figure('Renderer','painters','InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'FontSize',14);
plot0 = plot(t2_range_fs,to_plot_ground{3},'Parent',axes1,'LineWidth',2);
xlabel('t_2 (fs)','FontSize',16); ylabel('GS signal, a.u.','FontSize',14);
for lp=1:4
    set(plot0(lp),'DisplayName',strcat('[om_1,om_3]=[',num2str(om_pnts(lp,1)),...
                        ',',num2str(om_pnts(lp,2)),']cm^{-1}'));
end
set(plot0(lp),'Color',[0.749019622802734 0 0.749019622802734]);
legend1=legend(axes1,'show'); set(legend1,'FontSize',10);

%%
figure2 = figure('Renderer','painters','InvertHardcopy','off',...
    'Color',[1 1 1]);
axes2 = axes('Parent',figure2,'FontSize',14);
plot2 = plot(t2_range_fs,to_plot_excited{3},'Parent',axes2,'LineWidth',2);
xlabel('t_2 (fs)','FontSize',16); ylabel('ES signal, a.u.','FontSize',14);
for lp=1:4
    set(plot2(lp),'DisplayName',strcat('[om_1,om_3]=[',num2str(om_pnts(lp,1)),...
                        ',',num2str(om_pnts(lp,2)),']cm^{-1}'));
end
set(plot2(lp),'Color',[0.749019622802734 0 0.749019622802734]);
legend2=legend(axes2,'show'); set(legend2,'FontSize',10);

%%
for cmp = 1:5
    %%
om1_pnt = [1,1,2,2]; om3_pnt = [1,2,1,2]; %cmp =1
om_pnts=[12320,12320;12320,12520;12520,12320;12520,12520];
scale_fct2 = max(max(max(abs(R3(:,:,:,1)))));
to_plot1 = zeros(size(om_pnts,1),length(t2_range)); 
to_plot2 = to_plot1; to_plot3 = to_plot2;
for lp=1:4
    
    to_plot1(lp,:) = real(R3(om3_pnt(lp),:,om1_pnt(lp),cmp))/scale_fct2;
    to_plot2(lp,:) = real(R2(om3_pnt(lp),:,om1_pnt(lp),cmp))/scale_fct2;
    to_plot3(lp,:) = real(R5(om3_pnt(lp),:,om1_pnt(lp),cmp))/scale_fct2;                 
                 
end
%%
figure
plot0 = plot(t2_range_fs,to_plot1);
xlabel('t_2 (fs)'); ylabel('GS signal, a.u.');
for lp=1:4
    set(plot0(lp),'DisplayName',strcat('[om_1,om_3]=[',num2str(om_pnts(lp,1)),...
                        ',',num2str(om_pnts(lp,2)),']cm^{-1}'));
end
set(plot0(lp),'Color',[0.749019622802734 0 0.749019622802734]);
legend('show');
%%
figure
plot1 = plot(t2_range_fs,to_plot2+to_plot3);
xlabel('t_2 (fs)'); ylabel('ES signal, a.u.');
for lp=1:4
    set(plot1(lp),'DisplayName',strcat('[om_1,om_3]=[',num2str(om_pnts(lp,1)),...
                        ',',num2str(om_pnts(lp,2)),']cm^{-1}'));
end
set(plot1(lp),'Color',[0.749019622802734 0 0.749019622802734]);
legend('show');
end
%%
figure
plot(t2_range_fs,to_plot2)
xlabel('t_2 (fs)')
figure
plot(t2_range_fs,to_plot3)
xlabel('t_2 (fs)')
end