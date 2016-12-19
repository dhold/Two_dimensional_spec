%this should calculate the 2D spectra in the impulsive, RWA limit for the 
%Rephasing, non rephasing and coherent pathways.  
% S_k1(om_3,t_2,om_1) , S_k2(om_3,t_2,om_1) and S_k3(om3,om2,t1)

%% System related parameters
clear loaded_data
if ~exist('loaded_data','var') %skip this bit if data has already been loaded,
    %clear loaded_data in order to get this not to run
    for foldinglp =1:1 %pointless loop so I can fold code in editor
Temp = 300; %Kelvin
system_parameters_file = 'Parameters/dimer_system_params.mat';
[t_scale, B,speed_unit]= inv_cm_unit_sys(Temp);

%realisations of static disorder to generate
num_realisations = 1;

%% Beam params, wavevector and polarization
%specify a set of polarisations and wavevectors defining the beam
%propogation 

param_set = cat(3,[1/sqrt(2),0,1/sqrt(2);1/sqrt(2),0,1/sqrt(2);0,0,1;...
                 [0,1,0];[0,1,0];[1,1i,0]/sqrt(2);[1,-1i,0]/sqrt(2)],...
                 [1/sqrt(2),0,1/sqrt(2);1/sqrt(2),0,1/sqrt(2);0,0,1;...
                 [0,1,0];[0,1,0];[1,-1i,0]/sqrt(2);[1,+1i,0]/sqrt(2)]);
%param set contains in dim2 the spatial components of in dim 1
%k1,k2,k3,pol1,pol2,pol3,polhet.  dim3 contains all the other sets

%re/norephasing parameter ranges
%S(om_3,t_2,om_1) 
calc_rp = true;  calc_nr = true; 
om1_rng = linspace(10^7./(880),10^7./(700),100); %range 
t2_range = (0:5:1500)/1000*t_scale; 
om3_rng = linspace(10^7./(880),10^7./(700),100);

tauL = length(t2_range); w3L = length(om3_rng); w1L = length(om1_rng);

% Coherence related ranges (usually different)
%S(om_3,om_2,t_1) 
calc_coh = false; 

t1_coh = 40; %pass single to get the program to choose a range 
%based on markovian decay with this number of points (to exp(-5))
om2_coh= linspace(10^7./(880),10^7./(700),100); %range 
om3_coh= linspace(10^7./(880),10^7./(700),100); %range 

% times that will be used in the actual thing
t1_range = linspace(0,1/(om1_rng(2)-om1_rng(1)),...
   5*(om1_rng(end)-om1_rng(1))/(om1_rng(2)-om1_rng(1)));
t3_range = linspace(0,1/(om3_rng(2)-om3_rng(1)),...
    5*(om3_rng(end)-om3_rng(1))/(om3_rng(2)-om3_rng(1)));

%% Specify which signals to calculate

comp_to_calc = false(3,1);

%k_rp = k3+k2-k1;  %rephasing 
%k_nr = k3-k2+k1; %nom rephasing
if calc_rp || calc_nr
    comp_to_calc(1:2) = true;
end

% k_coh = -k3+k2+k1;%coherence
if calc_coh 
    comp_to_calc(3) = true;
end

%% System size and such

load(system_parameters_file,'H_site','sd_mat','BO_modes','Drude_modes',...
    'mu','R','mdip','mu_sd'); %loads parameters
loaded_data = true; 
N = length(H_site); %number of sites
mdip2 = +1i*pi*cross(R,mu); %other magnetic dipole factor, not include k factor

sd_noise = randn(num_realisations,length(sd_mat)); 
site_shift = sd_noise.*repmat(sd_mat,[num_realisations,1]); %sd_mat is diagonal uncorrelated noise
if ~isempty(mu_sd) %static disorder in mu,  assumed uncorrelated 
    %to site static disorder, stupid, but what can you do?
sd_noise2 = randn(num_realisations,length(sd_mat)); 
mu_shift = sd_noise2.*repmat(mu_sd,[num_realisations,1]);
end
%{
Drude_modes and BO_modes are size N cell arrays with each entry the
reorg E then damping /cutoff then om_0
lam_dru = Drude_modes{j}(:,1) , gam_dru = Drude_modes{j}(:,2)
lam_BO = BO_modes{j}(:,1) , gam_BO = BO_modes{j}(:,2), om0 ->(:,3)
 BO_modes{j}(:,4) = 0 by default, these can be changed in this section
 to include them in the Hamiltonian explicitly changing this number to a
 positive integer (which will remove them from the HEOM)
 can include any one of these modes on any site explicitly in the
 Hamiltonian 
 %}
explit_modes_in_ham = true; tot_vib = 1;

if explit_modes_in_ham 
    warning('not tested')
    %change the modes in the 4th coordinate to a positive integer so they
    %are explicitly included
    %BO_modes{1}(1,4) = 5; %would include a max of 5 excitations on this
    %mode for expample.
    BO_modes{1}(1,4) = 3;
    BO_modes{2}(1,4) = 3;
    lam_tot = zeros(N,1);
    for j = 1:N      
        inc_lg = BO_modes{j}(:,4)~=0;
        numvib{j} = BO_modes{j}(inc_lg,4);
        om_vib{j} =   BO_modes{j}(inc_lg,3);
        displ{j} = BO_modes{j}(inc_lg,1).^2./om_vib{j}; %check this is correct   
        tot_vib = tot_vib*prod(BO_modes{j}(:,4)); %total modes included
        lam_tot(j) = sum(BO_modes{j}(~inc_lg,1)) + sum(Drude_modes{j}(:,1));
    end    
    
    
else
    om_vib = []; numvib = []; displ = [];
    lam_tot = zeros(N,1);
    for k = 1:N
       lam_tot(k) = sum(BO_modes{k}(:,1)) + sum(Drude_modes{k}(:,1)) ;
    end  
    
end

[H_el,H_vib,H_e,H_f,fock_space_rep,M_p,M_pfull,~,~,H_coup] = ...
    generate_ex_vib_ham2(H_site,om_vib,numvib,displ,mu,mdip,R,[]); 

fock1 = fock_space_rep(2:N+1,:);
fock2 = fock_space_rep(N+2:end,:);

%{
% Generate the operator components of the interaction dipole moment 
%associated with each transition
V_n = zeros(1+N+N*(N-1)/2,1+N+N*(N-1)/2,N);
V_vec = zeros(1+N+N*(N-1)/2,1+N+N*(N-1)/2,3); %vector operator
V_mag = zeros(1+N+N*(N-1)/2,1+N+N*(N-1)/2,3); %higher order, mag dipole op
for lp =1:N
   
    V = zeros(1+N+N*(N-1)/2); 
    V(1,lp+1) = 1;  %mixes to ground
    
    lg = [false(N+1,1);fock_space_rep(N+2:end,lp)==1]; %elements with excitation at lp
    lg2 = fock_space_rep; lg2(:,lp) = 0; %[~,lg2] = find(lg2(lg,:)); 
    V(lg,2:N+1) =lg2(lg,:);
    V = V+V'; V_n(:,:,lp) = V;
    for j = 1:3
    V_vec(:,:,j) = V_vec(:,:,j) + mu(lp,j)*V_n(:,:,lp);
    V_mag(:,:,j) = V_mag(:,:,j) + m(lp,j)*V_n(:,:,lp);
    end
end

for j = 1:3
    V_ge(:,:,j) = kron(V_vec(:,:,j),eye(size(H_vib)));
end
%}

%% Precompute all the possible dipole expectation values 
% of the form <mu_4 dot e_out* ... mu_1 dot e_1  exp(i*k dot r_4..)>_{isotropic}
% these are computed as taylor expansions exp(1ix) ~ 1 + ix in terms of 4th
% and 5th order tensor averages

[av_4,av_5,av2_5]=ori_precomp_site_basis(param_set,mu,R,imag(mdip));

s_rf = [-1,+1,+1]; s_nr = [+1,-1,+1]; s_coh = [+1,+1,-1];
%these terms determine the order to sum in
    end
end

%% Next section deals with loops over static disorder, can be done parallel

calc_parallel = false;

if calc_parallel %pass parameters to paralell algorithm
    
    N_nodes = 10; 
    para_twoD_with_HEOM(N_nodes,num_realisations);
    
else 
for sdlp = 1:num_realisations

%%  Generate Hamiltonian things

om_dip = H_site + site_shift(sdlp);

%If this is slow I don't strictly need to generate all the vibrational
%stuff and coupling again, but currently this isn't a bottleneck
[H_el,H_vib,H_e,H_f,e_1,e_2,fock_space_rep,M_p,M_pfull,~,~,H_coup] = ...
    generate_ex_vib_ham2(om_dip,om_vib,numvib,displ,[],[],R,[]); 
M_e = M_p{1}; M_f = M_p{2};
%calculate the reorganisation energies in the exciton basis
lam_ex_full = zeros(1+N+N*(N-1)/2,1); 
 MM = blkdiag(1,M_e,M_f); 
for k = 1:N
    tmp = diag(mtimesx(MM,'C',mtimesx(diag(fock_space_rep(:,k)),MM)));
   lam_ex_full = lam_ex_full + tmp.*lam_tot(k);
end

%generate all the (vector) interaction operators in the exciton basis, both
%purely electric and electric + vib
[V_ge,V2_ge,V_ef,V2_ef,V_ge_ex,V2_ge_ex,V_ef_ex,V2_ef_ex...
,V_ge_full,V2_ge_full,V_ef_full,V2_ef_full,V_ge_ex_full,V2_ge_ex_full,V_ef_ex_full,V2_ef_ex_full]=...
  int_basis_ops('normal',fock_space_rep,mu,mdip+mdip2,M_p{1},M_p{2},tot_vib,M_pfull{1},M_pfull{2});
% specifying 'amp_angle' gives them in the form r, phi, theta which can be
% converted back via x  = r cos(phi) cos(theta), y  = r sin(phi)*cos(theta)
%
% alpha_n = mtimesx(MM,'C',mtimesx(V_n,MM));
% Cg  = squeeze(alpha_n(1,2:N+1,:)) ;  Cf = alpha_n(2:N+1,N+2:end,:);
% %coefficients mapping the interaction operator into the new basis
% 


%% Calculate exciton dipole moment averages and such like

%this is the "old" method, but not made a better one work yet...
mu_ex = squeeze(V_ge_ex); m_ex = squeeze(V2_ge_ex); 
mu_ex2 = permute(V_ef_ex,[2,3,1]); m_ex2 = permute(V2_ef_ex,[2,3,1]);
for j = 1:N
    m_ex(j,:) = m_ex(j,:)*e_1(j);
    for f = 1:N*(N-1)/2
        m_ex2(j,:,f) = m_ex2(j,:,f)*(e_2(f)-e_1(j));
    end
end


calc_OD  = true; %calculate off diagonal couplings
pol_set = param_set(4:7,:,1);
k_rp = param_set(3,:,1)+param_set(2,:,1) -param_set(1,:,1); %rephasing 
k_nr = param_set(3,:,1)-param_set(2,:,1) +param_set(1,:,1); %non rephasing

mag_rp = cross([param_set(1:3,:,1);k_rp],param_set(4:7,:,1));
mag_nr = cross([param_set(1:3,:,1);k_nr],param_set(4:7,:,1));

conj_rp = [true,false,false,true]; conj_nrp = [false,true,false,true];

[d_av_rp,d_av_rp_f,m_av_rp,m_av_rp_f]=...
        ori_av_fn(mu_ex,m_ex,mu_ex2,m_ex2,conj_rp,pol_set,mag_rp,calc_OD,false);

[d_av_nrp,d_av_nrp_f,m_av_nrp,m_av_nrp_f]=...
        ori_av_fn(mu_ex,m_ex,mu_ex2,m_ex2,conj_nrp,pol_set,mag_nr,calc_OD,false);    

%combine all the electric and magnetic dipoles into a combination of
%possible transition pathways that will eventually be considered
% [av_set1,av_set2,~] = ori_av_3rd_order2(mu_ex,m_ex,[],param_set,...
%                             [calc_rp,calc_nr,calc_coh],false);


%% Express transition operators in the site basis
V_ex_eg = zeros(N*tot_vib,tot_vib,N);
V_ex_fe = zeros(N*(N-1)/2*tot_vib,N*tot_vib,N,N*(N-1)/2);

temp_fe = zeros(N*(N-1)/2,N); temp_eg = zeros(N,1);
for j=1:N
    temp_eg2 = temp_eg; temp_eg2(j,1) = 1;
    V_ex_eg(:,:,j) = kron(M_e*temp_eg2,eye(tot_vib));
    for f=1:N*(N-1)/2
        
        temp_fe2 = temp_fe; temp_fe2(f,j) = 1;     
        temp_fe2 = M_f*temp_fe2*M_e'; 

        V_ex_fe(:,:,j,f) = kron(temp_fe2,eye(tot_vib));
    end
end

%% Calculate Redfield prop op and Lindblad op
[R_red,R_red_op] = redfield_calc_2(H_el(2:N+1,2:N+1),blkdiag(H_el(1,1),H_el(N+2:end,N+2:end)),...
    fock_space_rep,B,Drude_modes,[],BO_modes,[],[],true);
BO_modes{1}(1,4) = 3;   

[Lind_g,Lindblad_op_Lg] = Lindblad_op_gen(B,BO_modes,[],[],1,1);     
[Lind_e,Lindblad_op_Le] = Lindblad_op_gen(B,BO_modes,[],[],N,M_pfull{1});      
[Lind_f,Lindblad_op_Lf] = Lindblad_op_gen(B,BO_modes,[],[],N*(N-1)/2,M_pfull{2});      
        
[Lind_full,Lind_L] = Lindblad_op_gen(B,BO_modes,[],[],1+N+N*(N-1)/2,...
                        blkdiag(eye(tot_vib),M_pfull{1},M_pfull{2}));     
% *******************Next Calculate 2D signals******************

%% first calculate the initial coherences set up

%generate the N initial density matricies from the first pathways
prop_op_eg = H_prop_gen2(H_e,0,fock1,zeros(1,N),QQ_topass,const_factor,...
            coup_com_save,coup_acom_save);
om_scale_ge = e_1 + lam_ex_full(2:N+1);

t_to_pass = {t1_range,om1_rng,t1_range};
%t_to_pass = {t1_range,om1_rng,t1_coh}; %t1_range is range to calculate the whole function for and 
%then interpolate / fft for the other two variables
[rho_eg_coh,rho_eg_ft]=rho_eg_t_HEOM(V_ex_eg,1,N,tot_vib,t_to_pass,om_scale_ge, prop_op_eg,HL,HL);

%% Calculate the amplitude of the window functions

om_scale_ef = zeros(N,N*(N-1)/2);
for j=1:N
om_scale_ef(j,:) = E2 + lam_ex_full(N+2:end) - E1(j);
end

prop_op_fe = H_prop_gen2(H_f,H_e,fock2,fock1,QQ_topass,...
                          const_factor,coup_com_save,coup_acom_save);

om_scale_ge = E1 - lam_ex_full(2:N+1);

t_to_pass = {t3_range,-om3_rng};

[V_ge,V_ef] = window_fun_HEOM(V_ex_eg,V_ex_fe,N,H_trunc,HL,tot_vib,t_to_pass,...
                       prop_op_eg.',prop_op_fe.', -om_scale_ge,-om_scale_ef);
%reshape into Liouville space                                


%% Now calculate each of the required neq second order matricies

%om_scale_gf = E2 + lam_ex_full(N+2:end); %scaling for coherent cont only

% if calc_rp || calc_nr
% prop_op_ee = H_prop_gen2(H_e,H_e,fock1,fock1,QQ_topass,const_factor,coup_com_save,coup_acom_save);
% prop_op_gg = H_prop_gen2(0,0,zeros(1,N),zeros(1,N),QQ_topass,const_factor,coup_com_save,coup_acom_save);
% [rho_gg,rho_ee]=rho_neq_HEOM (V_ex_eg,N,1,H_trunc,HL,rho_eg_ft,prop_op_gg,prop_op_ee,t2_range) ;
% end
% if calc_coh
% prop_op_fg = H_prop_gen2(H_f,0,fock2,zeros(1,N),QQ_topass,const_factor,coup_com_save,coup_acom_save);
% simp_mix = true;
% [rho_fg]=rho_fg_HEOM (V_ex_fg,N,1,H_trunc,HL,rho_eg_coh,...
%                     sup_op_ge,t2_range,om2_rng,om_scale,simp_mix ) ;
% end

%this method works by calculating slices of the density matrix over t1 (or
%omega1 in frequency space) 
if calc_rp || calc_nr
prop_op_ee = H_prop_gen2(H_e,H_e,fock1,fock1,QQ_topass,const_factor,coup_com_save,coup_acom_save);
prop_op_gg = H_prop_gen2(0,0,zeros(1,N),zeros(1,N),QQ_topass,const_factor,coup_com_save,coup_acom_save);

sz_params = [N,HL,tot_vib]; %number of sites
[R1,R2,R3,R4,R5,R6] = rf_cont_HEOM(sz_params,rho_eg_ft,V_ex_eg,V_ex_fe,...
                            V_ge,V_ef,prop_op_gg,prop_op_ee,t2_range...
                            ,calc_rp,calc_nr,dip_av_set_rp,dip_av_set_nrp) ;
end
if calc_coh
    
om_scale_gf = E2 + lam_ex_full(N+2:end);    
    
prop_op_fg = H_prop_gen2(H_f,0,fock2,zeros(1,N),QQ_topass,const_factor,coup_com_save,coup_acom_save);
simp_mix = true;

[R7,R8]= coh_cont_HEOM (V_ex_fg,N,1,HL,V_ge,V_ef,rho_eg_coh,...
                    sup_op_ge,t2_range,om2_rng,om_scale,simp_mix ) ;
end

                                

%% Now calculate each of the different response function pathways


for j = 1:N  %first transition to exciton j and then either conjugate or j->f
    for j2 = 1:N %2nd
        for j3 = 1:N %3rd


        if calc_rp %235

    R2tmp = mtimesx(rho_ee(:,:,:,:,j,j2),'C',V_ex_eg(:,:,j3)); 
    %reshape into Liouville space to take inner product
    R2tmp = reshape(R2tmp,[N*H_trunc,size(rho_ee,4),size(rho_ee,5)]);
    
    R3tmp = mtimesx(V_ex_eg(:,:,j3),rho_gg(:,:,:,:,j,j2),'C'); 
    %reshape into Liouville space
    R3tmp = reshape(R3tmp,[N*H_trunc,size(rho_ee,4),size(rho_ee,5)]);
            for f= 1:N*(N-1)/2
    R5tmp = -mtimesx(V_ex_fe(:,:,j3,f),rho_ee(:,:,:,:,j,j2),'C'); 
    %reshape into Liouville space
    R5tmp = reshape(R5tmp,[N^2*(N-1)/2*H_trunc,size(rho_ee,4),size(rho_ee,5)]);    
            end
        end


        if calc_nr
    R1tmp = mtimesx(rho_ee(:,:,:,:,j),V_ex_eg(:,:,j2)); 
    %reshape into Liouville space
    R1tmp = reshape(R1tmp,[N*H_trunc_2,size(rho_ee,4),size(rho_ee,5)]);
    R1 = R1 + mtimesx(V_ge,'T',R1tmp); %Liouville space inner product
    
    R4tmp = mtimesx(V_ex_eg(:,:,j2),rho_gg(:,:,:,:,j)); 
    %reshape into Liouville space
    R4tmp = reshape(R4tmp,[N^2*(N-1)/2*H_trunc_2,size(rho_ee,4),size(rho_ee,5)]);
            for f= 1:N*(N-1)/2
    R6tmp = -mtimesx(V_ex_fe(:,:,j2),rho_ee(:,:,:,:,j)); 
    %reshape into Liouville space
    R6tmp = reshape(R6tmp,[N*tot_vib*H_trunc_2,size(rho_ee,4),size(rho_ee,5)]);    
            end
        end

        if calc_coh

    R7tmp = -mtimesx(rho_fg(:,:,:,:,j,f1),V_ex_eg(:,:,j2),'C'); 
    %reshape into Liouville space
    R7tmp = reshape(R7tmp,[N*tot_vib*H_trunc_2,size(rho_ee,4),size(rho_ee,5)]);
    
    R8tmp = mtimesx(V_ex_fe(:,:,j2),'C',rho_fg(:,:,:,:,j)); 
    %reshape into Liouville space
    R8tmp = reshape(R8tmp,[N^2*(N-1)/2*H_trunc_2,size(rho_ee,4),size(rho_ee,5)]);
    
        end
        end
    end
end


end

end