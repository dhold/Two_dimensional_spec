
N=2;
FWHM = [80,80,80]; %full width at half max of pulse intensity profile in fs
Peak_intensity = [0.5,0.5,0.5];
om_carrier = 8065.73*[2.21,2.11,2.11];
pol = [1,0,0;1,0,0;1,0,0];
beam_param = {FWHM,Peak_intensity,om_carrier,pol};
beam_param1 = {FWHM,[1,0,0].*Peak_intensity,om_carrier,pol};
beam_param2 = {FWHM,[0,1,0].*Peak_intensity,om_carrier,pol};
beam_param3 = {FWHM,[0,0,1].*Peak_intensity,om_carrier,pol};

%sig = FWHM/(4*sqrt(log(2)))*t_sc*1e-3; 
%E1env = @(t) Peak_intensity(1)*exp(-t.^2/sig(1)^2/2); 

load('Parameters/Hamiltonian_save.mat','PC645Hamiltonian','E0_PC645','PC645_dipole_data')

   sDO =  [0.5 0.1 0.05 0.04 0.04 0.02 0.03 0.01];
omega_DO = [10 200 665 818 1108 1270 1370 1645];

H_site = PC645Hamiltonian(1:N,1:N);
H_site = H_site - E0_PC645*(eye(N)); 
   
 %Set spectral density parameters 
 gam_dru = 100*ones(N,1); % drude decay constant / cut off frequency
 lam_dru = 100*ones(N,1); % renorm energy 
 
 clear B0_modes  Drude_modes 
 for lp=1:N;  Drude_modes{lp} = [lam_dru(lp),gam_dru(lp)]; end
 %Drude_modes{1} = zeros(0,2);  Drude_modes{2} = zeros(0,2); 
  gam_bro = 32; %damping
 mode_sel = [4]; %choose which modes to include in the spectral density
 nmodes = length(mode_sel); 
  %preassign empty
 for lp=1:N;  BO_modes{lp} = zeros(0,4); end
 
  for lp = mode_sel 
      BO_modes{1} = [BO_modes{1};[sDO(lp)*omega_DO(lp), gam_bro,omega_DO(lp),0]]; 
      BO_modes{2} = [BO_modes{2};[sDO(lp)*omega_DO(lp), gam_bro,omega_DO(lp),0]];
  end

 sys_params = {H_site,PC645_dipole_data,Drude_modes,BO_modes,[]};
 
 Kappa =1; max_tier = 2; use_RWA = true; om_s = 8065.73*2.11;
 tol = []; sd_mat = zeros(N); tau1 = 0; tau2 = 500;
 t_range_fs = linspace(-80*7,3500,2000); Temp = 300;
 [t_sc,B] = inv_cm_unit_sys(Temp);

%% Test linear polarization is as predicted

if 1==0
%very small peak intensity first beam only

[Pol_test,Pol_test2] = deal(zeros(3,length(t_range_fs)));

fle = open('Parameters/platonic_solid_rotations.mat');
%rot_set = fle.rot_set_iso_num;
rot_set = fle.rot_set_dodec_num;

x =[1,0,0]; y=[0,1,0]; z =[0,0,1]; 

pol_tmp = [x;x;x];  beam_param0 = {FWHM,[1e-4,0,0],om_carrier,pol_tmp}; 
last_rot = eye(3);
for lp = 1:length(rot_set)
    rot_mat = rot_set{lp}*last_rot; last_rot = rot_mat;
  Pol_test = Pol_test + non_perturbative_spec(Temp,beam_param0,sys_params,...
    t_range_fs, [0,0,0],tau1,tau2,rot_mat,max_tier,Kappa, use_RWA,om_s,sd_mat,tol);
  Pol_test2 = Pol_test2 + non_perturbative_spec(Temp,beam_param0,sys_params,...
    t_range_fs, [pi/2,0,0],tau1,tau2,rot_mat,max_tier,Kappa, use_RWA,om_s,sd_mat,tol);
end

%P^(1)(k,r) = P^+(t) e^(i k.r) + P^{-}(t) e^(-i k.r)
% so P^(1)(k,pi/2k) = P^+(t) i - P^{-}(t) i

figure; plot(t_range_fs, real(Pol_test - 1i*Pol_test2))
figure; plot(t_range_fs, real(Pol_test + 1i*Pol_test2)) 
figure; plot(t_range_fs, imag(Pol_test - 1i*Pol_test2))
end

%%

sig = (4*t_sc*1e-3)/(4*sqrt(log(2))); 

 t_range_fs_new = linspace(-80*7,1000,6500);
  t_range = t_range_fs_new/1e3*t_sc;
   dt = t_range(2)-t_range(1);  
  Om_rng = linspace(-1,1,length(t_range))*pi/dt;
  
  small_amp = 1e-4; %small amplitude
E1env = small_amp*exp(-t_range.^2/sig(1)^2/2); 

[~,lg] = min(abs(t_range_fs_new)); fft_disp = zeros(1,length(t_range_fs_new));
fft_disp(lg) = 1; fft_disp = fft(fft_disp);

[Pol_tmp,Pol_tmp2] = deal(zeros(3,length(t_range_fs_new),3));
x =[1,0,0]; y=[0,1,0]; z =[0,0,1];  %coordinate direction unit vectors
pol_tmp = [x;x;x];  rot_mat = eye(3); tol = [1e-20,1e-10];
clear rho_g rho_e rho_f
clever_prop = true; use_RWA = true;
tic
for lp = 1:3
    if lp == 2; pol_tmp = [y;x;x]; elseif lp==3; pol_tmp = [z;x;x]; end
    beam_param0 = {[400,eps,eps],[small_amp,0,0],om_carrier(1)*[1,1,1],pol_tmp};
  [Pol_tmp(:,:,lp)] = non_perturbative_spec(Temp,beam_param0,sys_params,...
    t_range_fs_new, [0,0,0],tau1,tau2,rot_mat,max_tier,Kappa, use_RWA,om_s,sd_mat,tol,clever_prop);
  Pol_tmp2(:,:,lp) =  non_perturbative_spec(Temp,beam_param0,sys_params,...
    t_range_fs_new, [pi/2,0,0],tau1,tau2,rot_mat,max_tier,Kappa, use_RWA,om_s,sd_mat,tol,clever_prop);

end
toc
%%

%scale out carrier frequency of the pulse in the forward propogating
%polarization
Pol_fwd = (Pol_tmp - 1i*Pol_tmp2).*repmat(exp(1i*om_carrier(1)*t_range),[3,1,3]);
Pol_fwd_iso = diagsum(Pol_fwd,1,3);
Pol_bc = (Pol_tmp + 1i*Pol_tmp2).*repmat(exp(-1i*om_carrier(1)*t_range),[3,1,3]);
Pol_bc_iso = diagsum(Pol_bc,1,3);
%take fourier transform

Pol_fwd_fft = fft(Pol_fwd_iso)./(3*fft_disp);  
Pol_fwd_fft = fftshift(Pol_fwd_fft); 

E1fft = fft(E1env)./fft_disp;   E1fft = fftshift(E1fft);

  figure; plot(Om_rng,real(Pol_fwd_fft)); hold on; plot(Om_rng,E1fft,'r')
  figure; plot(Om_rng,imag(Pol_fwd_fft)); hold on; plot(Om_rng,E1fft,'r')

%% undo the envelope
  

 % E1fft = fft( E1env.*cos((om_carrier(1)-om_mean)*t_range))./fft_disp; 

%  hold on; plot(Om_rng,E1fft * max(-real(S_iso)) / max(E1fft),'r');
  
  lg = abs(E1fft)>exp(-10)*max(E1fft); %reduced range that can be obtained from one pulse  
  
  %S(om-om_c) = P(om-om_c) / E(om)
  figure; plot(Om_rng(lg),real(Pol_fwd_fft(lg)./E1fft(lg)))
  %hold on; plot(Om_rng(lg),E1fft(lg)/max(E1fft),'r')
  figure; plot(Om_rng(lg),imag(Pol_fwd_fft(lg)./E1fft(lg)))
  %hold on; plot(Om_rng(lg),E1fft(lg)/max(E1fft),'r')  
 %% Method 1
 %choose phases so everything but the correct third order contributions
 %cancel.  
   %phases_set = {[0,0,0],[pi/2,0,0],[pi,0,0],[3*pi/2,0,0]};
   phases_set = {[0,0,0],[pi/2,0,0],[pi,0,0],[3*pi/2,0,0],...
                [0,0,pi],[pi/2,0,pi],[pi,0,pi],[3*pi/2,0,pi]};
   Pol_rp = zeros(3,length(t_range_fs)); Pol_nr = Pol_rp;
   cnst = [1,1i,-1,-1i]/4;   cnst2 = [1,-1i,-1,+1i]/4;  
   cnst = [cnst,-cnst];   cnst2 = [cnst2,-cnst2]; 
   
   phase_old = [inf,inf,inf]; %previous phases, choose impossible initial
 for lp = 1:length(phases_set)
     phases = phases_set{lp};
 [Pol_t] = non_perturbative_spec(Temp,beam_param,sys_params,...
    t_range_fs, phases,tau1,tau2,[],max_tier,Kappa, use_RWA,om_s,sd_mat,tol,true);
if phases(1) ~= phase_old(1) %don't both calculating if it is identical to before
  [Pol_1] = non_perturbative_spec(Temp,beam_param1,sys_params,...
    t_range_fs, phases,tau1,tau2,[],max_tier,Kappa, use_RWA,om_s,sd_mat,tol);
end
if phases(2) ~= phase_old(2)
  [Pol_2] = non_perturbative_spec(Temp,beam_param2,sys_params,...
    t_range_fs, phases,tau1,tau2,[],max_tier,Kappa, use_RWA,om_s,sd_mat,tol);
end
if phases(3) ~= phase_old(3)
  [Pol_3] = non_perturbative_spec(Temp,beam_param3,sys_params,...
    t_range_fs, phases,tau1,tau2,[],max_tier,Kappa, use_RWA,om_s,sd_mat,tol);
end

%Pol_save{lp} = {Pol_t,Pol_1,Pol_2,Pol_3};

Pol_rp = Pol_rp + cnst(lp)* (Pol_t - Pol_1 - Pol_2- Pol_3);
Pol_nr = Pol_nr + cnst2(lp)*(Pol_t - Pol_1 - Pol_2- Pol_3);

phase_old = phases;

 end
 
 %% Method 2, multidimensional DFT
 
 np = 4;  phase_range = linspace(0,2*pi,np+1);   

   Pol_rp = zeros(3,length(t_range_fs)); Pol_nr = Pol_rp;
   fct_rp = @(p1,p2,p3) exp(1i*(-p1+p2+p3));
   fct_nr = @(p1,p2,p3) exp(1i*(p1-p2+p3));
    fct_coh = @(p1,p2,p3) exp(1i*(p1+p2-p3));
    
 for lp1 = 1:np
     for lp2 = 1:np
         tic
         for lp3 = 1:np

     phases =  phase_range([lp1,lp2,lp3]);
 [Pol_t] = non_perturbative_spec(Temp,beam_param,sys_params,...
    t_range_fs, phases,tau1,tau2,[],max_tier,Kappa, use_RWA,om_s,sd_mat,tol);

Pol_rp = Pol_rp + fct_rp(phases(1),phases(2),phases(3))* (Pol_t);
Pol_nr = Pol_nr + fct_nr(phases(1),phases(2),phases(3))* (Pol_t);
         end
         toc
     end
 end
 
 %% Including orientation averaging
 
  x =[1,0,0]; y=[0,1,0]; z =[0,0,1];
required_pol = {[x;x;y;y],[y;y;x;x],[x;x;z;z],[z;z;x;x],[y;y;z;z],[z;z;y;y],...
           [x;y;x;y],[y;x;y;x],[x;z;x;z],[z;x;z;x],[y;z;y;z],[z;y;z;y],...     }
           [x;y;y;x],[y;x;x;y],[x;z;z;x],[z;x;x;z],[y;z;z;y],[z;y;y;z],...
            [x;x;x;x],[y;y;y;y],[z;z;z;z]}; %last line needed by all three

set_cont = {1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,[1,2,3],[1,2,3],[1,2,3]};

%pol_set = {x,y,z}; 
phi = (1+sqrt(5))/2; 
 v1 = [1,1,1]; v2 = [-1,1,1]; v3= [1,-1,1]; v4= [-1,-1,1];
v5 = [0,1/phi,phi]; v6 = [0,-1/phi,phi]; v7 = [1/phi,phi,0]; 
v8 = [-1/phi,phi,0]; v9 = [phi,0,1/phi]; v10 = [phi,0,-1/phi]; 
pol_set = {v1,v2,v3,v4,v5,v6,v7,v8,v9,v10}; 

num_pol = length(pol_set);
[Pol_1 , Pol_2, Pol_3] = deal(zeros(3,length(t_range_fs),num_pol));

for lp = 1:num_pol  %calculate the one pulse averages
    pol_lp = pol_set{lp}; pol_lp = pol_lp/norm(pol_lp);
    tmp = beam_param1; tmp{4} = pol_lp;
    Pol_1(:,:,lp) = non_perturbative_spec(Temp,tmp,sys_params,...
    t_range_fs, [0,0,0],tau1,tau2,[],max_tier,Kappa, use_RWA,om_s,sd_mat,tol,true);
    tmp = beam_param2; tmp{4} = pol_lp;
    Pol_2(:,:,lp) = non_perturbative_spec(Temp,tmp,sys_params,...
    t_range_fs, [0,0,0],tau1,tau2,[],max_tier,Kappa, use_RWA,om_s,sd_mat,tol,true);
    tmp = beam_param3; tmp{4} = pol_lp;
    Pol_3(:,:,lp) = non_perturbative_spec(Temp,tmp,sys_params,...
    t_range_fs, [0,0,0],tau1,tau2,[],max_tier,Kappa, use_RWA,om_s,sd_mat,tol,true);
end

[Pol_1av , Pol_2av, Pol_3av] = deal(zeros(1,length(t_range_fs)));

for lp = 1:num_pol  %calculate the one pulse averages
    pol_lp = pol_set{lp}; pol_lp = pol_lp/norm(pol_lp);
 Pol_1av = Pol_1av + Pol_1(1,:,lp).*pol_lp(1)+Pol_1(2,:,lp).*pol_lp(2)+Pol_1(3,:,lp).*pol_lp(3);
 Pol_2av = Pol_2av + Pol_2(1,:,lp).*pol_lp(1)+Pol_2(2,:,lp).*pol_lp(2)+Pol_2(3,:,lp).*pol_lp(3);
 Pol_3av = Pol_3av + Pol_3(1,:,lp).*pol_lp(1)+Pol_3(2,:,lp).*pol_lp(2)+Pol_3(3,:,lp).*pol_lp(3);
end
 Pol_1av = Pol_1av / num_pol;   Pol_2av = Pol_2av / num_pol; 
  Pol_3av = Pol_3av / num_pol; 
  
  
  %% Method 4, choose molecular rotations that should give average
  
     phases_set = {[0,0,0],[pi/2,0,0],[pi,0,0],[3*pi/2,0,0],...
                [0,0,pi],[pi/2,0,pi],[pi,0,pi],[3*pi/2,0,pi]};
   % phases_set = {[0,0,0],[pi/2,0,0],[pi,0,0],[3*pi/2,0,0]};
   cnst = [1,1i,-1,-1i]/4;   cnst2 = [1,-1i,-1,+1i]/4;  
   cnst = [cnst,-cnst];   cnst2 = [cnst2,-cnst2]; 
  

phi = (1+sqrt(5))/2;  %golden ratio
%list verticies of a regular icosahedron
%vert_pos = [0,1,phi;0,-1,phi;0,1,-phi;0,-1,-phi]; %cyclic perms of these
%vert_pos = [vert_pos;vert_pos(:,[3,1,2]);vert_pos(:,[2,3,1])]/sqrt(1+phi^2); 
%list verticies of a regular dodecehedron   
vert_pos = [1,1,1;-1,1,1;1,-1,1;-1,-1,1];
vert_pos = [vert_pos;-vert_pos];
v_tmp = [0,1/phi,phi; 0,-1/phi,phi ; 0,1/phi,-phi ; 0,-1/phi,-phi];
vert_pos = [vert_pos;v_tmp;v_tmp(:,[3,1,2]);v_tmp(:,[2,3,1])];
%this is a set of 20 vectors pointing out from the centre of a dodecehdron
%through to the verticies

Pol_rp = zeros(3,length(t_range_fs),length(vert_pos)); Pol_nr = Pol_rp;

  x =[1,0,0]; y=[0,1,0]; z =[0,0,1]; pols = [x;x;x];

rot_set = {eye(3)}; %first rotation is the identity
%solve V2 = R*V for all different R -> V2*V^(-1) = R
X = linsolve(vert_pos(1,:).',eye(3)); %gives V^(-1)
for lp = 2:length(vert_pos)
    rot_set{lp} = vert_pos(lp,:).' * X;
end
    %pol_lp = [x;x;x]; %set at the start anyway
    clever_prop = false; %prop individual components seperately when 
    %electric field is weak, saves some time when I have tested it
    
for lp = 1:length(rot_set)  %calculate the one pulse averages
    tic
    rot_mat = rot_set{lp};

    Pol_1(:,:,lp) = non_perturbative_spec(Temp,beam_param1,sys_params, t_range_fs,...
    [0,0,0],tau1,tau2,rot_mat,max_tier,Kappa, use_RWA,om_s,sd_mat,tol,clever_prop );

    Pol_2(:,:,lp) = non_perturbative_spec(Temp,beam_param2,sys_params,t_range_fs, ...
    [0,0,0],tau1,tau2,rot_mat,max_tier,Kappa, use_RWA,om_s,sd_mat,tol,clever_prop );
init_phase_3 = tau2*om_carrier(3);
    Pol_3(:,:,lp) = non_perturbative_spec(Temp,beam_param3,sys_params,t_range_fs, ...
     [0,0,init_phase_3],tau1,tau2,rot_mat,max_tier,Kappa, use_RWA,om_s,sd_mat,tol,clever_prop );

 %{
%     for lp2 = 1:length(phases_set) %this calculates the three pulse averages
%         phase_lp = phases_set{lp2}; cc = exp(1i*phase_lp); 
%         Pol_full = non_perturbative_spec(Temp,beam_param,sys_params,t_range_fs,...
%          phase_lp,tau1,tau2,rot_mat,max_tier,Kappa, use_RWA,om_s,sd_mat,tol,clever_prop );        
%         Pol_rp(:,:,lp) = Pol_rp(:,:,lp) + cnst(lp2)* (Pol_full-...
%         cc(1)*Pol_1(:,:,lp)+cc(2)*Pol_2(:,:,lp)+cc(3)*Pol_3(:,:,lp));
%         Pol_nr(:,:,lp) = Pol_nr(:,:,lp) + cnst2(lp2)*(Pol_full-...
%         cc(1)*Pol_1(:,:,lp)+cc(2)*Pol_2(:,:,lp)+cc(3)*Pol_3(:,:,lp));
%     end
   %}
     one_loop_time = toc;
end

%Pol_rp_av = sum(Pol_rp,3); Pol_nr_av = sum(Pol_nr,3);
% for lp2 = 1:length(phases_set) %average and remove the individual components
%     phase_lp = phases_set{lp2}; cc = exp(1i*phase_lp); %phase of indiv pulses
%     Pol_rp_av = Pol_rp_av - cnst(lp2)*sum(cc(1)*Pol_1+cc(2)*Pol_2+cc(3)*Pol_3,3);
%     Pol_nr_av = Pol_nr_av - cnst(lp2)*sum(cc(1)*Pol_1+cc(2)*Pol_2+cc(3)*Pol_3,3);
% end



%% Hetdetect

FWHM_het = 160; %full width at half max of pulse intensity profile in fs
om_het_rng = 8065.73*(2.01+linspace(-0.05,0.05));
pol = [1,0,0];


 %%
 pol_to_plot = sum(Pol_1,3);%Pol_rp_av;
  figure; plot(t_range_fs, pol_to_plot); 
  
  lg = t_range_fs >0;
 T_range = t_range_fs(lg)/1e3*t_sc; 
  om_HET = 8065.73*2.01; %range of heterodyne frequencies
 % scale_fc = repmat(exp(1i*t_range*om_HET),[3,1]);
 %figure; plot(t_range_fs, imag(scale_fc.*pol_to_plot)); 
 
   dt = T_range(2)-T_range(1); 
  %T = T_range(end)-T_range(1); 
  om_rng = linspace(-1,1,length(T_range))*pi/dt;
 %scale phase oscillations
 Pol_ft = fftshift(fft(pol_to_plot(:,lg),[],2),2);

 figure; plot(om_rng,real(Pol_ft))
 
 %%
 
 [tmp1,tmp2,tmp3] = deal(t_range_fs*0);  [~,lg] = min(abs(t_range_fs)); 
 [~,lg2] = min(abs(t_range_fs -tau1));   [~,lg3] = min(abs(t_range_fs -tau1-tau2)); 
 tmp1(lg) = 1; tmp2(lg2) = 1; tmp3(lg3) = 1;
 Tmp1 = fft(tmp1); Tmp2 = fft(tmp2); Tmp3 = fft(tmp3);
 
 sig = FWHM/(4*sqrt(log(2)))*t_sc*1e-3; 
E1env = Peak_intensity(1)*exp(-t_range.^2/sig(1)^2/2); 
 
  t_range = t_range_fs/1e3*t_sc;
   dt = t_range(2)-t_range(1);    T = t_range(end)-t_range(1); 
  om_rng2 = linspace(-1,1,length(t_range))*pi/dt;

  
 scale_fc = repmat(exp(-1i*t_range*(om_s)),[3,1]);
 
  pol_to_plot = Pol_1.*scale_fc;%Pol_rp.*scale_fc;%Pol_rp_av;
 figure; plot(om_rng2, fftshift(fft(sum(pol_to_plot)/3)./Tmp1))
 %plot(om_rng2, [fftshift(fft(pol_to_plot,[],2)./repmat(Tmp1,3,1),2)]);%...
                            %;0*fftshift(fft(E1env)./Tmp1)])
  
  %%
 freq_range = linspace(-2000,4000,5000);
test = zeros(3,length(freq_range)); 
 
env_fn = 1;%exp(-((t_range_fs-500).^2/500^2/2));

 for lp = 1:length(freq_range)
 scale_fc = repmat(env_fn.*exp(1i*t_range*freq_range(lp)),[3,1]);
 test(:,lp) = freq_range(lp)*trapz(t_range,scale_fc.*pol_to_plot,2);
 end
 
 figure; plot(freq_range,imag(test));
 
 %include envelope
 env_fn = exp(-freq_range.^2/100^2/2);
 c = conv_fft(imag(test(1,:)), env_fn,'same');
  figure; plot(freq_range,c)