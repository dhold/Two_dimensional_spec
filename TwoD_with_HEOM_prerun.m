function [coup_com_save,coup_acom_save,const_factor,QQ_topass,nn,fock_space_rep,...
    H_site,site_shift,g_cnst,lam_tot,V_n,av_2,av_3,av2_3,av_4,av_5,av2_5,chk] = ...
    TwoD_with_HEOM_prerun(system_parameters_file,Temperature_in_Kelvin,...
    num_realisations,max_tier,Kappa,param_set) %#ok<*STOUT>

%This prerun function does the calculations which are independent of static
%disorder. This includes the calculation of the coupling between tiers of
%the HEOM, averages of the form <mu_1 mu_2 m_3 mu_4> 

%calculate scaling from ps->inv cm, thermo beta in cm^(-1) 
[t_scale, B,speed_unit]= inv_cm_unit_sys(Temperature_in_Kelvin);


%% load system parameters from file
load(system_parameters_file,'H_site','sd_mat','BO_modes','Drude_modes',...
    'mu','R','mdip','mu_sd'); %loads parameters
N = length(H_site); %number of sites
if all(mu_sd==0) %#ok<*NODEF>
    mu_sd = []; end

sd_noise = randn(num_realisations,length(sd_mat)); 
site_shift = sd_noise.*repmat(sd_mat,[num_realisations,1]); %sd_mat is diagonal uncorrelated noise

if ~isempty(mu_sd) %static disorder in mu,  assumed uncorrelated 
    %to site static disorder
    warning('currently static disorder in mu is not supported')
    sd_noise2 = randn(num_realisations,length(sd_mat)); 
    mu_shift = sd_noise2.*repmat(mu_sd,[num_realisations,1]);
end
%% Calculate HEOM components coupling tiers
clear cc vv cc2 cc_com cc_acom

g_cnst = zeros(N,1); lam_tot = zeros(N,1); %For Markov approx
    for j = 1:N
        if ~isempty(BO_modes{j}) %check if underdamped modes present
lam_tot(j) = sum(BO_modes{j}(:,1)) + sum(Drude_modes{j}(:,1)) ;
        %include modes with last entry = 0
        to_inc = BO_modes{j}(:,4)==0;
    om_0 =   BO_modes{j}(to_inc,3); %#ok<*USENS>
    gam  =   BO_modes{j}(to_inc,2);
    lambda = BO_modes{j}(to_inc,1);
        else
lambda = []; gam = []; om_0 = [];   %assume all empty
lam_tot(j)=sum(Drude_modes{j}(:,1)) ; 
        end
if size(Drude_modes{j},2)==1 %this allows the passing of just a pure dephasing term
    gam_dru = []; lam_dru =[];
    QQ_extra = [Drude_modes{j}(1),Drude_modes{j}(2)]; 
    %note passing a second argument will include dissipation 
else
    gam_dru = Drude_modes{j}(:,1); lam_dru=  Drude_modes{j}(:,2);
    QQ_extra = [0,0];
end    
    
%calculate markovian linebroadening function to estimate decay rates
    g_cnst(j) = line_broad_fn_markov(B,gam,om_0,lambda,gam_dru,lam_dru);

%This function calculates the coefficients for the commutator and anti
%commutator coefficients in the HEOM, along with the truncation correction
    [cc1,cc2R,cc2I,vv1,vv2,QQ,cc_tmp] = coeffients_from_brownian_new...
   (lambda,gam,om_0,Temperature_in_Kelvin,Kappa,gam_dru,lam_dru); 
    if ~isempty(cc_tmp)
     cc_com{j}= [cc2R,cc1]; vv{j} = [vv2,vv1]; cc(j,:)=cc_tmp; %#ok<*AGROW>
     cc_acom{j}= [cc2I,cc1*0]; else
      cc_com{j}= []; vv{j} = [];cc_acom{j}= [];
    end
    QQ_topass(j,:) = [QQ,0]+QQ_extra;
    end
[coup_com_save,coup_acom_save,const_factor,nn]=...
    multi_site_drude_and_brownian_HOM_4(cc_com,cc_acom,vv,inf,max_tier);



%% Make a matrix relating which site the excited state is present on 
% for each state in the hilbert space

fock1 = eye(N);
fock_space_rep = zeros(1+N+N*(N-1)/2,N); fock_space_rep(2:N+1,:) = fock1;
cnt = 0;
for j = 1:N
        rng = j+1:N; 
    for k = 1:length(rng)
        cnt = cnt+1;  kk = rng(k);
       fock_space_rep(N+1+cnt,j) = 1; fock_space_rep(N+1+cnt,kk) = 1; 
    end
end

% Generate the operator components of the interaction dipole moment 
%associated with each transition
V_n = zeros(1+N+N*(N-1)/2,1+N+N*(N-1)/2,N);
for lp =1:N
   
    V = zeros(1+N+N*(N-1)/2); 
    V(1,lp+1) = 1;  %mixes to ground
    
    lg = [false(N+1,1);fock_space_rep(N+2:end,lp)==1]; %elements with excitation at lp
    lg2 = fock_space_rep; lg2(:,lp) = 0; %[~,lg2] = find(lg2(lg,:)); 
    V(lg,2:N+1) =lg2(lg,:);
    V = V+V'; V_n(:,:,lp) = V;
    
end

%% Precompute all the possible dipole expectation values in site basis
% of the form <mu_4 dot e_out* ... mu_1 dot e_1  exp(i*k dot r_4..)>_{isotropic}
% these are computed as taylor expansions exp(1ix) ~ 1 + ix in terms of 4th
% and 5th order tensor averages

if ~isempty(param_set)

if ~iscell(param_set)
[~,~,~,av_4,av_5,av2_5]=ori_precomp_site_basis(param_set,mu,R,imag(mdip));
chk =[];
else %compute different order with different param_set
[~,~,~,av_4,av_5,av2_5,chk]=ori_precomp_site_basis2(param_set{1},param_set{2},mu,R,imag(mdip));
end
%also compute some first order averages with CP light
 tmp = [[1,+1i,0]/sqrt(2);[1,-1i,0]/sqrt(2);0,0,1];
 
[av_2,av_3,av2_3]=ori_precomp_site_basis(tmp,mu,R,imag(mdip));
 
end
