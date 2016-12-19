for jj =1:3
clear params_for_wrapper
T = 77;  [~ , B ,~ ]= inv_cm_unit_sys(T); %needed to convert parameters in 
%their equation into parameters for the HEOM

E_a = 12328;
E_b = 12472;
J = 70.7;
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
 lam_bro =  lam_bro*sqrt(2)
 
 gam_bro = gamma_mode;
 BO_modes{1} = [lam_bro, gam_bro,om_0_rng(jj),0];
 BO_modes{2} = [lam_bro, gam_bro,om_0_rng(jj),0];

gamma_site = 0.025*delta_Ex; %pure dephasing for each site
%relate the pure dephasing for each site to a Drude spectral density
%2 \lambda / \beta \gamma, assume \gamma = 100 and thus
%  \lambda = \beta \gamma gamma_site / 2
if 1 ==0
gam_dru = 100;
lam_dru = gam_dru*B*gamma_site/2; %gives an unconvincingly low number :/
%lam_dru = lam_dru*2*pi %increase by factor of 2*pi
Drude_modes{1} = [lam_dru,gam_dru];
Drude_modes{2} = [lam_dru,gam_dru];
else %just include some pure dephasing shit
    %passing in this format will be pure dephasing
   % Drude_modes{1} = [gamma_site*2*pi;0]; 
   % Drude_modes{2} = [gamma_site*2*pi;0];
    Drude_modes{1} = [4*gamma_site;0]; 
    Drude_modes{2} = [4*gamma_site;0];
end

mu = [1, 0.5, 0;0, 1, 0];
mu_ex = M_prj*mu;
mu_sd = [0,0]; sd_mat = [0,0]; %no static disorder parameters
R = [1, 0, 0;0, 0, 0]; mdip = mu*0; %beyond dipole shit, ignored
pdm_shift = []; %no shift to double excited state
%% Save
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
    save(strcat('saved_data/plenio_w_pd_paramset_',num2str(jj),'.mat'),...
        'R1','R2','R3','R4','R5','R6','lin_spec','tmp_gg_t','tmp_ee_t','numvib',...
    'Drude_modes','BO_modes','opts');    
    
end
%%
%figure
 for jj =1:3 %plot the graphs from plenios paper
    
    load(strcat('saved_data/plenio_w_pd_paramset_',num2str(jj),'.mat'),'R3')

scale_fct2 = max(max(max(abs(real(R3(:,:,:,1))))));
pnt = length(t2_range);
%subplot(3,1,jj)
figure
pcolor(om3_rng,om1_rng,abs(real(squeeze(R3(:,pnt,:,1)))).'/scale_fct2)
shading interp
xlabel('\omega_3 (cm^{-1})')
ylabel('\omega_1 (cm^{-1})')
 end
%%
 [a,pnt2] = min(abs(t2_range - 1.6/delta_Ex)); %find min 
load(strcat('saved_data/plenio_w_pd_paramset_',num2str(2),'.mat'),'R3','R2','R5')
 %%
%pnt2=pnt2+1 %iterate by 1 to see tdep, pretty cool!

%load(strcat('test_data.mat'),'R3','R2','R5')
figure
pcolor(om3_rng,om1_rng,squeeze(abs(real(R3(:,pnt2,:,1)+R2(:,pnt2,:,1)+R5(:,pnt2,:,1)))).'/scale_fct2)
shading interp
xlabel('\omega_3')
ylabel('\omega_1')

%% also plot the slices
om_pnts=[12320,12320;12320,12520;12520,12320;12520,12520]; %om1,om3
%find closest points to these
to_plot1 = zeros(size(om_pnts,1),length(t2_range)); 
to_plot2 = to_plot1; to_plot3 = to_plot2;
for lp=1:4
    [~,om1_pnt]=min(abs(om_pnts(lp,1)-om1_rng));
    [~,om3_pnt]=min(abs(om_pnts(lp,2)-om3_rng));
    to_plot1(lp,:) = real(R3(om3_pnt,:,om1_pnt,1))/scale_fct2;
    to_plot2(lp,:) = real(R2(om3_pnt,:,om1_pnt,1))/scale_fct2;
    to_plot3(lp,:) = real(R5(om3_pnt,:,om1_pnt,1))/scale_fct2;                 
                 
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
if 1==0
%% Next bit

params_for_wrapper = 'Parameters/file_for_wrapper.mat';
Temperature_in_Kelvin = 77; system_parameters_file= 'parameters/plenio_paper_parameters.mat';
t1_coh =[]; om2_coh=[]; om3_coh = []; %no coherence shit
pol = [1,0,0]; kvec = [0,0,1];
beam_param_set = [kvec;kvec;kvec;pol;pol;pol;pol];
%same shit as before
E_a = 12328;E_b = 12472;J = 70.7;H_site = [E_a,J;J,E_b];
[M_prj,E_ex]= eig(H_site);E_ex = diag(E_ex);  delta_Ex = E_ex(2)-E_ex(1);
gamma_mode  = 0.005*delta_Ex; SHR = 0.02; om_0_rng = delta_Ex; 
  om_0 = {om_0_rng,om_0_rng}; lam_bro = SHR*om_0_rng; 
 gam_bro = gamma_mode; BO_modes{1} = [lam_bro, gam_bro,om_0_rng,0];
 BO_modes{2} = [lam_bro, gam_bro,om_0_rng,0];
gamma_site = 0.025*delta_Ex; 
Drude_modes = {[gamma_site*2*pi;0],[gamma_site*2*pi;0]}; ;
mu = [1, 0.5, 0;0, 1, 0];mu_ex = M_prj*mu;    mu_sd = [0,0]; 
R = [1, 0, 0;0, 0, 0]; mdip = mu*0; pdm_shift = [];
sd_mat = [0.17*delta_Ex,0.17*delta_Ex];  %DOES contain static D
%much smaller range of frequencies to consider
om1_rng = [12320,12520]; om3_rng = om1_rng; t2_range_fs = linspace(0,2000,210);
numvib = {3,3};
num_realisations = 50;

save(params_for_wrapper,'Temperature_in_Kelvin','system_parameters_file',...
 'beam_param_set','om1_rng','t2_range_fs','om3_rng','t1_coh','om2_coh','om3_coh','num_realisations')
%%
%twoD_with_HEOM_wrapper;
%save('saved_data/plenio_sd_inc2.mat','R1','R2','R3','R4','R5','R6'...
%        ,'lin_spec','tmp_gg_t','tmp_ee_t','max_tier','Kappa','num_realisations');
  twoD_with_Redfield_wrapper;
  save('saved_data/plenio_w_exp_mode2.mat','R1','R2','R3','R4','R5','R6'...
          ,'lin_spec','tmp_gg_t','tmp_ee_t','numvib','num_realisations');

%%
om1_pnt = [1,1,2,2]; om3_pnt = [1,2,1,2];
om_pnts=[12320,12320;12320,12520;12520,12320;12520,12520];
scale_fct2 = max(max(max(abs(R3(:,:,:,1)))));
to_plot1 = zeros(size(om_pnts,1),length(t2_range)); 
to_plot2 = to_plot1; to_plot3 = to_plot2;
for lp=1:4
    
    to_plot1(lp,:) = real(R3(om3_pnt(lp),:,om1_pnt(lp),1))/scale_fct2;
    to_plot2(lp,:) = real(R2(om3_pnt(lp),:,om1_pnt(lp),1))/scale_fct2;
    to_plot3(lp,:) = real(R5(om3_pnt(lp),:,om1_pnt(lp),1))/scale_fct2;                 
                 
end
%%
figure
plot(t2_range_fs,to_plot1)
xlabel('t_2 (fs)')
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
%%
figure
plot(t2_range_fs,to_plot2)
xlabel('t_2 (fs)')
figure
plot(t2_range_fs,to_plot3)
xlabel('t_2 (fs)')
end