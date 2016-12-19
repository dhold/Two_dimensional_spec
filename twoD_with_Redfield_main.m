function [R1_sdav,R2_sdav,R3_sdav,R4_sdav,R5_sdav,R6_sdav,lin_spec,tmp_gg_t,tmp_ee_t]=...
        twoD_with_Redfield_main(sdlp_rng,fock_space_rep,H_el,H_vib,H_e,H_f,...
       site_shift,B,Drude_modes,BO_modes,V_n,V2_n,av_2,av_3,av2_3,av_4,av_5,av2_5,  ...
        calc_rp,calc_nr,om1_rng,t2_range,om3_rng,opts,t_verbose);

    
%this should calculate the 2D spectra in the impulsive, RWA limit for the 
%Rephasing, non rephasing and coherent pathways.  
% S_k1(om_3,t_2,om_1) , S_k2(om_3,t_2,om_1) and S_k3(om3,om2,t1)
%%
 
if isempty(av_4) && isempty(av_5) && isempty(av2_5)
    mu =av_2; compute_dipole_moments_in_ex_basis = true
    R = av_3{1}; mm = av_3{2};
    param_set = av2_3;
else
    compute_dipole_moments_in_ex_basis = false;
end
    
if ~exist('opts','var')
   opts = 'Redfield'; 
elseif ~any([strcmp(opts,'Redfield'),strcmp(opts,'Forster'),strcmp(opts,'Pure_dephasing')])
    error('opts should be one of either Redfield, Forster or Pure_dephasing')
end

N = size(fock_space_rep,2);  N2 = N*(N-1)/2; 
Nv = length(H_vib); sz_full = (1+N+N2)*Nv;

fock1 = fock_space_rep(2:N+1,:);
fock2 = fock_space_rep(N+2:end,:);
sdlp=1;
if ~exist('t_verbose','var') %give value true to output time taken 
    t_verbose = false;
end
%tauL = length(t2_range); w3L = length(om3_rng); w1L = length(om1_rng);
%%
for sdlp = sdlp_rng

%%  Generate Hamiltonian things
% add static disorder
H_site_shift1 = site_shift(sdlp,:)*(fock1.'); %shift to ham
H_ee = H_el(2:N+1,2:N+1) + diag(H_site_shift1); [P_e,E1] = eig(H_ee); 
H_site_shift1 = kron(diag(H_site_shift1),eye(size(H_vib))); %pad out to full size
H_e = H_e + H_site_shift1; 

H_site_shift2 = diag(site_shift(sdlp,:)*(fock2.')); %same for double excited
H_ff = H_el(N+2:end,N+2:end) + H_site_shift2; [P_f,E2] = eig(H_ff);  
H_site_shift2 = kron(H_site_shift2,eye(size(H_vib)));
H_f = H_f + H_site_shift2;

 E1= diag(E1);  E2= diag(E2); 

%P_e_pad = kron(P_e,eye(size(H_vib)));
%P_f_pad = kron(P_f,eye(size(H_vib)));

MM = blkdiag(1,P_e,P_f);
MM_full = kron(blkdiag(1,P_e,P_f),eye(size(H_vib)));
%projectors in site only basis for single and excited states

%projection interaction operator from exciton basis ->dipole basis
alpha_n = mtimesx(MM,mtimesx(V_n,MM,'C'));

Cg  = squeeze(alpha_n(1,2:N+1,:)) ;  Cf = alpha_n(2:N+1,N+2:end,:);
%coefficients mapping the interaction operator into the new basis

%% Calculate Dephashing tensor
if strcmp(opts,'Redfield') 

sz1 =  1+length(H_e) + length(H_f);  Nv =  length(H_vib);

R_red = redfield_calc_2(H_e,blkdiag(1,H_f),...
            fock_space_rep,B,Drude_modes,[],BO_modes,[],[],true);      
%modes with "zero" in the 4th column will be included in the rates in the
%Redfield tensor
        
LL = sz1^2*Nv^2;
R_red_op = zeros(LL); 

temp = zeros(LL);
cnt1 = 1; cnt2 = 1;
cnt1vib =0; cnt2vib =1; 
%warning('this bit seems wrong')
%cnt_test = 0;
for lp = 1:LL^2
    
    %iterate along one vibrational lvl
    cnt1vib = cnt1vib+1;
    if cnt1vib > Nv  %if beyond this level reset to one and add another to cnt1
        cnt1vib=1; cnt1 = cnt1+1;
        if cnt1 > sz1
            cnt1 = 1;  cnt2vib = cnt2vib+1;
            if cnt2vib > Nv
                cnt2vib = 1; cnt2 = cnt2+1;
            end
        end
    end
    temp = temp*0;
   for a = 1:sz1 %loop over electronic degrees of freedom
        for b = 1:sz1 
            temp((a-1)*Nv+cnt1vib,(b-1)*Nv+cnt2vib) = R_red(cnt1,cnt2,a,b);
        end
   end
%cnt_test = cnt_test + sum(sum(temp~=0));
    R_red_op(lp,:) = reshape(temp,1,sz1^2*Nv^2 ).';

end
    R_red_op = sparse(R_red_op);
    
elseif strcmp(opts,'Forster') 
    error('Forster not included yet')
    
elseif strcmp(opts,'Pure_dephasing') 
        %includes a [V,[V,rho]] style term, still given the same variable
        %name
        
        R_red_op = sparse([],[],[],sz_full^2,sz_full^2);
        
      for j = 1:length(Drude_modes)
          
          QQ = Drude_modes{j}; %Drude_modes will be two numbers, one with a
          %pure dephasing part and another with a dissipative coefficient
          %site basis interaction operator to the exciton basis
             V = sparse(MM_full*kron(diag(fock_space_rep(:,j)),speye(Nv))*MM_full');
           Qcom =  kron(speye(size(V)),V)-kron(V.',speye(size(V)));
           Qacom = kron(speye(size(V)),V)+kron(V.',speye(size(V))) ; %anti commutator
           
            R_red_op = R_red_op + QQ(1).*Qcom*Qcom + QQ(2).*Qacom*Qcom;
          
      end      
      R_red_op = sparse(R_red_op);
end



%% Calculate dipole moments averages in the exciton basis and such like


if ~compute_dipole_moments_in_ex_basis

%1st order averages
av_set_fo = mtimesx(mtimesx(Cg,'C',av_2),Cg); 
av_set_cd = mtimesx(mtimesx(Cg,'C',av_3+av2_3),Cg); 

%3rd order averages

tmp = mtimesx(mtimesx(Cg,'C',av_4),Cg); %first two transitions are between
%g and one exciton transitions
tmp = permute(tmp,[3,4,1,2,5]); %permute so next two dimensions are in mtimesx

%the first one involves transitions between g-e only (GSB, SE)
av_set_GSB = mtimesx(mtimesx(Cg,'C',squeeze(tmp)),Cg);
av_set_GSB = permute(av_set_GSB,[3,4,1,2,5]); %permute back to original order
%the next loop is more complicated because it involves transitions from e-f
%and then f'-e', need to reshape
Cf_flat = reshape(Cf,N^2*(N-1)/2,N); %flatten that shit out

av_set_ESA = mtimesx(mtimesx(Cf_flat,'C',squeeze(tmp)),Cf_flat); %ESA type
av_set_ESA = permute(av_set_ESA,[3,4,1,2,5]);
%av_set_ESA (e1,e2,e3*(N-1)*N/2+f1,e4*(N-1)*N/2+f2) type calls
%next calculate the higher order moments from the other interactions

tmp = mtimesx(mtimesx(Cg,'C',av_5+av2_5),Cg);
tmp = permute(tmp,[3,4,1,2,5,6]); 

%approximate amplitude of k1, k2 by the vertical transition energy of the
%exciton states this map into k~w_res/c
%next I want to really scale SE and GSB differently to account for Stokes
%shift, this is a closer approximation Do this for each path when
%performing the loop
% E_GSB = E1+lam_ex_full(2:N+1);  
% E_SE = E1-lam_ex_full(2:N+1);  

av_set2_GSB  = mtimesx(mtimesx(Cg,'C',tmp), Cg);
av_set2_GSB = permute(av_set2_GSB,[3,4,1,2,5,6]);
av_set2_ESA  = mtimesx(mtimesx(Cf_flat,'C',tmp), Cf_flat);  
av_set2_ESA = permute(av_set2_ESA,[3,4,1,2,5,6]);

else

[C,~] = eig(H_el(2:N+1,2:N+1));
[CC,~] = eig(H_el(N+2:end,N+2:end));
fock2 = fock_space_rep(N+2:end,:);

mu_e_f = zeros(N*(N-1)/2,3,N); mu_ex = zeros(N,3);
for k = 1:N
    for n = 1:N
    mu_ex(k,:) = mu_ex(k,:) + C(n,k)*mu(n,:); %this is what I think is
    %correct as I need the participation of the nth site in the kth exciton
    %the wavefunction for the kth exciton is 
    %|psi_k> = sum_n C(n,k) |psi_n>
    %mu_ex(k,:) = mu_ex(k,:) + C(k,n)*mu(n,:);
    end
    for f = 1:N*(N-1)/2

for n=1:N
    for m = n+1:N
        lg = logical(fock2(:,n))&logical(fock2(:,m));
        mu_e_f(f,:,k) =  mu_e_f(f,:,k) + CC(lg,f)*(C(n,k)*mu(m,:) + C(m,k)*mu(n,:));
       % mu_e_f(f,:,k) =  mu_e_f(f,:,k) + CC(f,lg)*(C(k,n)*mu(m,:) + C(k,m)*mu(n,:)); 
    end
end
    end
end

av_set_GSB = zeros(N,N,N,N,size(param_set,3));
av_set_ESA = zeros(N,N,N^2*(N-1)/2,N^2*(N-1)/2,size(param_set,3));
av_set2_GSB = zeros(N,N,N,N,size(param_set,3),3);
av_set2_ESA  = zeros(N,N,N^2*(N-1)/2,N^2*(N-1)/2,size(param_set,3),3);
av_set_cd = zeros(N,N,3);  av_set_fo = zeros(N,N); 
for l = 1:size(param_set,3)
    pol = param_set(4:7,:,l);
for k1 = 1:N %slow loop
    for k2 = 1:N
        
 av_set_fo(k1,k2) = dot(mu_ex(k1,:),mu_ex(k2,:))/3; %lowest order
        
       for k3 = 1:N 
           for k4 = 1:N  
              
              mu_set = [mu_ex(k1,:);mu_ex(k2,:);mu_ex(k3,:);mu_ex(k4,:)];
              av_set_GSB(k1,k2,k3,k4,l) = tensor_av(mu_set, pol);
              for f=1:N*(N-1)/2
                    
                  mu_set = [mu_ex(k1,:);mu_ex(k2,:);mu_e_f(f,:,k3);mu_e_f(f,:,k4)];
                  av_set_ESA(k1,k2,(k3-1)*N*(N-1)/2+f,(k4-1)*N*(N-1)/2+f,l)...
                            = tensor_av(mu_set, pol);

              end
           end
       end
    end
end
end

end
% Next Calculate 2D signals

%% calculate operator coupling excitations on each site to the exciton
%hamiltonian
                
H_tot = blkdiag(H_vib,H_e,H_f); H_ex =  MM_full'*H_tot*MM_full;                          
L_full =  -1i*(kron(eye(length(H_ex)),H_ex)-kron(H_ex.',eye(length(H_ex))));
Lindblad_op = Lindblad_op_gen(B,BO_modes,[],[],1+N+N*(N-1)/2,MM_full);

prop_op_full = L_full + Lindblad_op - R_red_op;

    tmpg = zeros(sz_full); tmpg(1:Nv,1:Nv)=1;
	tmpg = logical(reshape(tmpg,sz_full^2,1)); %picks this out

%reduced operators acting only on ground excited coherences, these don't
%mix to p_gg or p_ee'

    tmpge = zeros(sz_full); tmpge(1:Nv,Nv+1:Nv*(N+1))=1;
    tmpge = logical(reshape(tmpge,sz_full^2,1));
    tmpeg = zeros(sz_full); tmpeg(Nv+1:Nv*(N+1),1:Nv)=1;
    tmpeg = logical(reshape(tmpeg,sz_full^2,1));
    
    tmpe  = zeros(sz_full); 
    tmpe(Nv+1:Nv*(N+1),Nv+1:Nv*(N+1)) = 1;
	tmpe = logical(reshape(tmpe,sz_full^2,1));
    
%Select particular sections of this

prop_op_eg = prop_op_full(tmpeg,tmpeg);
prop_op_ge = prop_op_full(tmpge,tmpge);
prop_op_ee = prop_op_full(tmpe,tmpe);
prop_op_gg = prop_op_full(tmpg,tmpg);
%% first calculate the initial coherences set up

rho_eq = diag(exp(-B*diag(H_vib))); rho_eq = rho_eq/trace(rho_eq);
rho_eq = reshape(rho_eq,[numel(rho_eq),1]);

%generate the N initial density matricies from the first pathways            
                
 tic
 [rho_eg_ft,rho_ge_ft]=rho_ge_f_Redfield(N,Nv,rho_eq,om1_rng,prop_op_eg,prop_op_ge);
%rho_eg(om_1) and rho_ge(-om_1)
init_coherence_time_to_calc = toc;
if t_verbose
    init_coherence_time_to_calc
end
%prop_op_freq = (1i*om1_rng*eye(size(prop_op_eg))-prop_op_eg)^(-1);

 %dt =  t1_range(2)-t1_range(1); om1 = pi*(-1/dt:2/t1_range(end):1/dt);

%use this coherence value to calculate linear spectra 

sigma_om1 = zeros(1,size(rho_eg_ft,2)); alpha_om1 = zeros(1,size(rho_eg_ft,2));
for j=1:N
    for j2 = 1:N
            V_sec = zeros(Nv,N*Nv); V_sec(1:Nv,((j-1)*Nv+1):j*Nv) = eye(Nv);
            V_sec = reshape(V_sec.',[numel(V_sec),1]).';
        sigma_om1 = sigma_om1 + av_set_fo(j,j2)*V_sec*rho_eg_ft(:,:,j); %lowest order
        alpha_om1 = alpha_om1 + av_set_cd(j,j2)*V_sec*rho_eg_ft(:,:,j); %cd / OR
    end
end
alpha_om1 = alpha_om1.*reshape(om1_rng,size(alpha_om1)); %also include scaling
lin_spec{1,:}={sigma_om1,alpha_om1};
%% Calculate the amplitude of the window functions

%{
om_scale_ef = zeros(N,N*(N-1)/2);
for j=1:N
om_scale_ef(j,:) = E2 - E1(j);% + lam_ex_full(N+2:end)
end
V_coup2 = zeros(N/2*(N-1),N/2*(N-1),N); %system bath coupling
for j = 1:N %project each term into the exciton basis
V_coup2(:,:,j) = M_f'*diag(fock2(:,j))*M_f;
end
 for j = 1:N %I think this should be the same but am not 100% sure
     coup_com_save2{j} = coup_com_save{j}.'; 
     coup_acom_save2{j} = coup_acom_save{j}.'; 
 end   
 
 %when H_prop_gen2 is given a 9th argument it assumes the time dependence
 %is to the left.  One passes the same for |f><e| operator time dep as you
 %would for |e><f| density matrix time dep and the function switches
prop_op_eg2 = H_prop_gen2(diag(E1),0,V_coup,zeros(1,1,N),QQ_topass,...
                         const_factor,coup_com_save2,coup_acom_save2,1);
prop_op_fe2 = H_prop_gen2(diag(E2),diag(E1),V_coup2,V_coup,QQ_topass,...
                          const_factor,coup_com_save2,coup_acom_save2,1);
%}
tic 
[V_ge,V_ef]=window_fun_Redfield_f(N,Nv,om3_rng,prop_op_full);
time_dep_of_output_op_calc_time = toc;
if t_verbose
    time_dep_of_output_op_calc_time
end
%This is expressed in Liouville space                                

%also use this time dependent value to calculate linear spectra 
% dt3 =  t3_range(2)-t3_range(1); om3 = pi*(-1/dt3:2/t3_range(end):1/dt3);

sigma_om3 = zeros(size(V_ge,1),1); alpha_om3 = zeros(size(V_ge,1),1);
for j=1:N
    V_sec = zeros(N*Nv,Nv); V_sec((j-1)*Nv+1:j*Nv,1:Nv) = eye(Nv);
    V_sec = kron(speye(Nv),V_sec);
    rho_sec = V_sec*rho_eq;
    
    for j2 = 1:N
        sigma_om3 = sigma_om3 + av_set_fo(j,j2)*V_ge(:,:,j2)*rho_sec;
        alpha_om3 = alpha_om3 + av_set_cd(j,j2)*V_ge(:,:,j2)*rho_sec;
    end
end
alpha_om3 = alpha_om3.*reshape(om3_rng,size(alpha_om3));
lin_spec{2,:}={sigma_om3,alpha_om3};


%% Now calculate each of the required neq second order matricies

%this method works by calculating slices of the density matrix over t1 (or
%omega1 in frequency space) 
if calc_rp || calc_nr

sz_params = [N,Nv];

[R1,R2,R3,R4,R5,R6,tmp_gg_t,tmp_ee_t] = rp_cont_Redfield2(sz_params,rho_eg_ft,...
                    rho_ge_ft,V_ge,V_ef,prop_op_gg,prop_op_ee,t2_range,E1,E2...
       ,calc_rp,calc_nr,av_set_GSB,av_set_ESA,av_set2_GSB,av_set2_ESA,[],t_verbose) ;   
   
   
else
    R1 = []; R2 =[]; R3 = []; R4 = []; R5 = []; R6 = [];
end


%% Add together to create average

if sdlp == 1
    
    R1_sdav = R1;  R4_sdav = R4;   R6_sdav = R6;   
    R2_sdav = R2;  R3_sdav = R3;   R5_sdav = R5; 
else %add extra terms
    R1_sdav = R1_sdav + R1;  R4_sdav = R4_sdav + R4;   R6_sdav = R6_sdav + R6;   
    R2_sdav = R2_sdav + R2;  R3_sdav = R3_sdav + R3;   R5_sdav = R5_sdav + R6; 
end
end