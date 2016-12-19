%% Initialise parameters

params_for_wrapper = 'Parameters/file_for_wrapper.mat';
Temperature_in_Kelvin = 300; 
system_parameters_file= 'parameters/dimer_CD_parameters.mat';
file_name = 'saved_data/10degset_HEOM_pt40pi.mat'; %file name to save under

[t_scale, B, ~]= inv_cm_unit_sys(Temperature_in_Kelvin); 

%derive expressions for the polarizations
clear beam_param_set
syms a b c a2 b2 c2
ry = [cos(b),0,sin(b);0,1,0;-sin(b),0,cos(b)];
rz = [cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
rz2 = [cos(c),-sin(c),0;sin(c),cos(c),0;0,0,1];
 
rot_op = rz2*ry*rz;
rot_op2 = subs(rot_op,[a,b,c],[a2,b2,c2]);

init_rot =pi/4*[1,3,1,-3];
%undo rotation as my code assumes the het detect is along z
 
b_sub = 10*pi/180;

 fb_l = @(b) b.^2/4+b.^4/24 + (0.6840*b.^6  -1.7634*b.^8-0.6270*b.^10)*1e-3; 
 fb_u = @(b) 3*b.^2/4+15*b.^4/24 + 0.4146*b.^6  +  0.2183*b.^8  -0.3670*b.^10;            
 
 cnt = 0; bpset1_save = zeros(4,3,length(b_sub)*2);
        bpset2_save = zeros(7,3,length(b_sub)*2);
 for lp = 1:length(b_sub)  
     for lp2 = 1:2
cnt = cnt + 1; %count parameter
t1 = init_rot(1); t2 = init_rot(2); t3 = init_rot(3); t4 = init_rot(4);
if lp2 == 2
t1=t1-fb_u(b); t2 = t2 + fb_u(b); t3 = t3 - fb_l(b); t4 = t4 + fb_l(b);   
end
e3 = subs(rot_op*[1;0;0],[c,a],[pi/4,t3]);
e4 = subs(rot_op*[1;0;0],[c,a],[3*pi/4,t4]);
e1 = subs(rot_op*[1;0;0],[c,a],[7*pi/4,t1]);
e2 = subs(rot_op*[1;0;0],[c,a],[5*pi/4,t2]);

k1 = subs(rot_op*[0;0;1],[c,a],[7*pi/4,t1]); %+kx - ky
k2 = subs(rot_op*[0;0;1],[c,a],[5*pi/4,t2]); %-kx - ky
k3 = subs(rot_op*[0;0;1],[c,a],[pi/4,t3]);   %+kx + ky
k4 = subs(rot_op*[0;0;1],[c,a],[3*pi/4,t4]); %-kx + ky


%In this code het detection is assumed polarized along x and
%with wavevector z, rotate coordinates as required
%%
%rotate to new basis as my code assumes the het detect is along z polarized
%along x, also make row vectors because again the code assumes this.
rot_tmp = subs(rot_op2,[a2,b2,c2],[-3*pi/4 , -b,-t4]);
rot_tmp = subs(rot_tmp,b,b_sub);

k1=(double(rot_tmp)*double(subs(k1,b,b_sub))).'; e1=(double(rot_tmp)*double(subs(e1,b,b_sub))).';
k2=(double(rot_tmp)*double(subs(k2,b,b_sub))).'; e2=(double(rot_tmp)*double(subs(e2,b,b_sub))).';
k3=(double(rot_tmp)*double(subs(k3,b,b_sub))).'; e3=(double(rot_tmp)*double(subs(e3,b,b_sub))).';
k4=(double(rot_tmp)*double(subs(k4,b,b_sub))).'; e4=(double(rot_tmp)*double(subs(e4,b,b_sub))).';

bpset1_save(:,:,cnt) = [e1;e2;e3;e4];
bpset2_save(:,:,cnt) = [e1;e2;e3;e4;k1;k2;k3];
     end
 end
%beam configurations which will be required
 x=[1,0,0]; y=[0,1,0]; z=[0,0,1]; 

beam_param_set{1} = bpset1_save;%[e1;e2;e3;e4]; % coherence selec
beam_param_set{2} = bpset2_save;%[e1;e2;e3;e4;k1;k2;k3];
  F4 = ori_F_calc(beam_param_set{1});
  F5= ori_F_calc(cat(3,[e1;e2;e3;e4;x],[e1;e2;e3;e4;y],[e1;e2;e3;e4;z]));

%parameters for the system 

%exciton parameters
%%
%theta_ex_rng = linspace(0,1/2,1000)*pi;
%for k = 1:1000
    
E_L = 12000; E_U = 12400; DE = E_U-E_L;  theta_ex = 0.40*pi;%theta_ex_rng(k);
%derive site parameters from these exciton parameters
J = DE/sqrt(4+(tan(theta_ex)-1/tan(theta_ex))^2) %J_save(k) = J; 

if mod(theta_ex/pi,1/2) < 0.25;  J = -J; end
deltaE = sqrt(DE^2-4*J^2); %DE_save(k) = deltaE;
E_a = (E_L+E_U - deltaE)/2; E_b = (E_L+E_U + deltaE)/2;
%E_a = 12000;    E_b = 12400;   J = 150;   

H_site = [E_a,J;J,E_b];

[M_prj,E_ex]= eig(H_site);     E_ex = diag(E_ex); 
%these two should be the same as the originals
delta_Ex = E_ex(2)-E_ex(1); theta_ex2 = atan(-M_prj(1,1)./M_prj(2,1));
M_prj_test = sign(J)*[-sin(theta_ex2),cos(theta_ex2);cos(theta_ex2),sin(theta_ex2)];
%end
%%
mu = [1, 0.5, 0;0, 1, 0];   mu_ex = M_prj*mu;    mu_sd = [0,0]; 
R = 10^(-7)*[0, 1, 0;  0, 0, 1]; mdip = mu*0; pdm_shift = [];
sd_mat = [0,0];  %No static D

%gamma_site = 0.025*delta_Ex;  %site based pure dephasing rate
l_dr = 30; g_dru = 100;
%%
%range of frequencies to consider

om1_rng = linspace(E_ex(1)-300,E_ex(2)+300,250);% %range in cm^(-1) of {freq,time}
om3_rng = linspace(E_ex(1)-300,E_ex(2)+300,360);
t2_range_fs = 0;  %just consider zero delay spectra

 num_realisations = 1;
 if num_realisations == 1; sd_mat = sd_mat*0; end %no static disorder

   % Drude_modes = {[sqrt(l_dr*g_dru);0],[sqrt(l_dr*g_dru);0]}; %pure dehphasing 
    Drude_modes = {[l_dr,g_dru],[l_dr,g_dru]}; 
    BO_modes{1} =[]; BO_modes{2}= [];

save(system_parameters_file,'R','mu','H_site','BO_modes',...
        'Drude_modes','sd_mat','mdip','mu_sd','pdm_shift')


t1_coh = []; om2_coh = []; om3_coh = [];  pop_only = false; numvib = 1;
max_tier = 4; Kappa = 2; 
save(params_for_wrapper,'Temperature_in_Kelvin','system_parameters_file','numvib',...
 'Kappa','max_tier','pop_only','beam_param_set','om1_rng','t2_range_fs','om3_rng',...
 't1_coh','om2_coh','om3_coh','num_realisations')    
    

%%
    twoD_with_HEOM_wrapper;
%save in a format that can be loaded easily
save(file_name,'R1','R2','R3','R4','R5','R6','om1_rng','om3_rng','t2_range_fs','Drude_modes','BO_modes'...
        ,'lin_spec','tmp_gg_t','tmp_ee_t','max_tier','Kappa','num_realisations','beam_param_set','chk',...
        'b_sub','H_site','-v7.3');


%%    
    
CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
CMRmap=interp1(1:9,CMRmap,linspace(1,9,18));
%xlabel('\tau (fs)','FontSize',16); ylabel('Probe frequency \omega_3 cm^{-1}','FontSize',14);

NR_sig_im = squeeze(imag(R1+R4+R5)); RF_sig_im = squeeze(imag(R2+R3+R6));
NR_sig_re = squeeze(real(R1+R4+R5)); RF_sig_re = squeeze(real(R2+R3+R6));

%sf = max(max(abs(NR_sig_im(:))),max(abs(RF_sig_im(:))));
sf = 1e-07
%correct for code being... well wrong for the chiral part
len = size(NR_sig_im,3)/2;
NR_sig_im(:,:,len+1:end) = NR_sig_im(:,:,len+1:end)*2*pi; 
RF_sig_im(:,:,len+1:end) = RF_sig_im(:,:,len+1:end)*2*pi;

tp = 2;
if 1 == 0  %change 0 to 1 to plot all these
figure; pcolor(om1_rng,om3_rng,abs(NR_sig_im(:,:,tp)/sf));  colormap(CMRmap); shading flat
xlabel('\omega_1/(2 \pi c) cm^{-1}'); ylabel('\omega_s/(2 \pi c) cm^{-1}');
figure; pcolor(om1_rng,om3_rng,abs(RF_sig_im(:,:,tp)/sf));  colormap(CMRmap); shading flat
xlabel('\omega_1/(2 \pi c) cm^{-1}'); ylabel('\omega_s/(2 \pi c) cm^{-1}');

figure; pcolor(om1_rng,om3_rng,abs(NR_sig_im(:,:,tp+len)/sf));  colormap(CMRmap); shading flat
xlabel('\omega_1/(2 \pi c) cm^{-1}'); ylabel('\omega_s/(2 \pi c) cm^{-1}');
figure; pcolor(om1_rng,om3_rng,abs(RF_sig_im(:,:,tp+len)/sf));  colormap(CMRmap); shading flat
xlabel('\omega_1/(2 \pi c) cm^{-1}'); ylabel('\omega_s/(2 \pi c) cm^{-1}');

figure; pcolor(om1_rng,om3_rng,abs(sum(NR_sig_im(:,:,[tp,tp+len]),3)/sf));  colormap(CMRmap); shading flat
xlabel('\omega_1/(2 \pi c) cm^{-1}'); ylabel('\omega_s/(2 \pi c) cm^{-1}');
figure; pcolor(om1_rng,om3_rng,abs(sum(RF_sig_im(:,:,[tp,tp+len]),3)/sf));  colormap(CMRmap); shading flat
xlabel('\omega_1/(2 \pi c) cm^{-1}'); ylabel('\omega_s/(2 \pi c) cm^{-1}');
end
%%
%figure; pcolor(om1_rng,om3_rng,sum(RF_sig_im,3)/sf);  colormap(CMRmap); shading flat
%xlabel('\omega_1/(2 \pi c) cm^{-1}'); ylabel('\omega_s/(2 \pi c) cm^{-1}');

figure1 = figure('Renderer','painters',...
    'PaperSize',[20.98404194812 29.67743169791],...
    'InvertHardcopy','off', 'Colormap',CMRmap, 'Color',[1 1 1]);

axes1 = axes('Parent',figure1,'FontSize',16);

box(axes1,'on');   hold(axes1,'all');
% Create surface
contourf('Parent',axes1,om1_rng,om3_rng,sum(RF_sig_im(:,:,[tp,tp+len]),3)/sf,length(CMRmap));

% Create xlabel
xlabel('\omega_1/(2 \pi c) cm^{-1}','FontSize',16);

% Create ylabel
ylabel('\omega_s/(2 \pi c) cm^{-1}','FontSize',16);

% Create colorbar
colorbar('peer',axes1,'FontSize',16);

if 1 == 1
%%
for k = 1:2
      %  if k == 1; fnm ='10deg_pt38pi.mat';  fnm2 ='10deg_pt39pi.mat';
      %  elseif k == 2;  fnm ='10deg_nocomp_pt38pi.mat';
      %           fnm2 ='10deg_nocomp_pt39pi.mat';
      %  end
      fnm ='10degset_HEOM_pt38pi';  fnm2 ='10degset_HEOM_pt40pi';
load(fnm,'R2','R3','R6')
 RF_sig_im = squeeze(imag(R2+R3+R6)); len = size(RF_sig_im,3)/2;
load(fnm2,'R2','R3','R6')
 RF_sig_im = RF_sig_im - squeeze(imag(R2+R3+R6)); 
 
sf = 1e-07;  RF_sig_im(:,:,k+len) = RF_sig_im(:,:,k+len)*2*pi;
figure1 = figure('PaperSize',[20.98404194812 29.67743169791],...
    'InvertHardcopy','off', 'Colormap',CMRmap, 'Color',[1 1 1]); 
axes1 = axes('Parent',figure1,'FontSize',16);
box(axes1,'on');   hold(axes1,'all');
contourf('Parent',axes1,om1_rng,om3_rng,sum(RF_sig_im(:,:,[k,k+len]),3)/sf,length(CMRmap));
xlabel('\omega_1/(2 \pi c) cm^{-1}','FontSize',16);
ylabel('\omega_s/(2 \pi c) cm^{-1}','FontSize',16);

end
end

if 1 == 0
%%

      fnm ='10degset_HEOM_pt38pi';  fnm2 ='10degset_HEOM_pt40pi';
load(fnm,'R2','R3','R6')
 RF_sig_im = squeeze(imag(R2+R3+R6)); len = size(RF_sig_im,3)/2;
load(fnm2,'R2','R3','R6')
 RF_sig_im2 = squeeze(imag(R2+R3+R6)); 
 
sf = 1e-07;  RF_sig_im(:,:,k+len) = RF_sig_im(:,:,k+len)*2*pi;
RF_sig_im2(:,:,k+len) = RF_sig_im2(:,:,k+len)*2*pi;

toplot1= sum(RF_sig_im(:,:,[1,1+len]),3)/sf; toplot2 = sum(RF_sig_im2(:,:,[1,1+len]),3)/sf;
toplot3= sum(RF_sig_im(:,:,[2,2+len]),3)/sf; toplot4 = sum(RF_sig_im2(:,:,[2,2+len]),3)/sf;

min_c_1 = min(min(toplot1(:)),min(toplot2(:))); max_c_1 = max(max(toplot1(:)),max(toplot2(:)));
min_c_2 = min(min(toplot3(:)),min(toplot4(:))); max_c_2 = max(max(toplot3(:)),max(toplot4(:)));

figure; 
subplot(2,2,1); contourf(om1_rng/10^4,om3_rng/10^4,toplot1,length(CMRmap)); caxis([min_c_1, max_c_1])
subplot(2,2,2); contourf(om1_rng/10^4,om3_rng/10^4,toplot2,length(CMRmap)); caxis([min_c_1, max_c_1])
subplot(2,2,3); contourf(om1_rng/10^4,om3_rng/10^4,toplot3,length(CMRmap)); caxis([min_c_2, max_c_2])
xlabel('\omega_1/(2 \pi c)  [10^{4} cm^{-1}]','FontSize',16);
ylabel('\omega_s/(2 \pi c)  [10^{4} cm^{-1}]','FontSize',16);
subplot(2,2,4); contourf(om1_rng/10^4,om3_rng/10^4,toplot4,length(CMRmap)); caxis([min_c_2, max_c_2])
%colormap(CMRmap)

figure; 
subplot(1,2,1); contourf(om1_rng/10^4,om3_rng/10^4,toplot1-toplot2,length(CMRmap)); 
subplot(1,2,2); contourf(om1_rng/10^4,om3_rng/10^4,toplot3-toplot4,length(CMRmap));
end