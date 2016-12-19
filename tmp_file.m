
T = 77;  [t_sc, B,~]= inv_cm_unit_sys(T); %needed to convert parameters in 
%their equation into parameters for the HEOM

E_a = 10000;
E_b  = E_a ; J = -300;
H_site = [E_a,J;J,E_b];

[M_prj,E_ex]= eig(H_site);
E_ex = diag(E_ex);  delta_Ex = E_ex(2)-E_ex(1);
%delta_Ex = delta_Ex*2*pi % is this where they went wrong / I did??
      
 BO_modes = {[],[]};
gam_dru = 10/t_sc; %1/100fs in inv cm
lam_dru = 60;
Drude_modes = {[lam_dru,gam_dru],[lam_dru,gam_dru]};

mu = [5, 0, 0;-1, 0, 0];
mu_ex = M_prj*mu;
mu_sd = [0,0]; sd_mat = [0,0]; %no static disorder parameters
R = [1, 0, 0;-1, 0, 0]; mdip = mu*0; %beyond dipole shit, ignored

pdm_shift{1,2} = 200; pdm_shift{2,1} = 200;
%% Save
save('Adv_HEOM_for_Efficient_eval_params.mat','R','mu','H_site','BO_modes',...
        'Drude_modes','sd_mat','mdip','mu_sd','pdm_shift')