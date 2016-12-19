function [GSB_sav,SE_sav,ESA_sav,GSB_cav,SE_cav,ESA_cav,rho_gg_t,rho_ee_t,lin_spec]=...
        twoD_main_new(fock_space_rep,H_el,H_vib,H_e,H_f,...
       site_shift,B,Drude_modes,BO_modes,QQ,const_factor,V_n,V_n2,av_2,av_3,av2_3,av_4,av_5,av2_5,  ...
        t1_range,om1_rng,t2_range,om3_rng,opts,t_verbose);
    
%This calculates the 2D signal S_R/NR(om3,t2_rng,om1 or t1)
% the probe detection is assumed to be frequency
% resolved, so it is sufficient to calculate individual frequencies without
% knowledge of the pulse shape.  
% This function uses either the HEOM or Redfield / pure dephasing for the 
% propogators (depending on opts). If tau = 0 or [] then the code works in
% frequency space for the initial coherence and final interaction operator,
% accounting for the frequency shifts due to coherent oscillation in the 
% population time.  It works in either the site, exciton or exciton
% vibrational basis.
% Here "Drude_modes" and "BO_modes" are only passed as Drude_modes and
% BO_modes if Redfield is used, otherwise Drude/BO_modes should be passed as
% Coup_com_save and Coup_acom_save. QQ is the pure dephasing term used in 
% the HEOM, note that the pure dephasing is just the HEOM with 0 tiers.
% It also either works in the exciton basis or exciton vibrational basis
% depending on the options of 
%%
N = size(fock_space_rep,2);
fock1 = fock_space_rep(2:N+1,:);
fock2 = fock_space_rep(N+2:end,:);
if strcmp(opts(2),'HEOM')
HL = size(const_factor,1); %length of HEOM
coup_com_save = Drude_modes; coup_acom_save = BO_modes;
clear Drude_modes BO_modes
else
    HL = 1; 
end
sdlp=1;
if ~exist('t_verbose','var') %give value true to output time taken 
    t_verbose = false;
end
%tauL = length(t2_range); w3L = length(om3_rng); w1L = length(om1_rng);
%%
for sdlp = sdlp_rng

%%  Diagonalise Hamiltonian blocks for projection operators

%include static D
H_site_shift1 = site_shift(sdlp,:)*(fock1.'); %shift to ham
H_ee = H_el(2:N+1,2:N+1) + diag(H_site_shift1); 
H_site_shift1 = kron(diag(H_site_shift1),eye(size(H_vib))); %pad out to full size
H_e = H_e + H_site_shift1; 

H_site_shift2 = diag(site_shift(sdlp,:)*(fock2.')); %same for double excited
H_ff = H_el(N+2:end,N+2:end) + H_site_shift2; 
H_site_shift2 = kron(H_site_shift2,eye(size(H_vib)));
H_f = H_f + H_site_shift2;

 %diagonalise the blocks of the electronic Hamiltonian
[P_e,e1] = eig(H_ee);  [P_f,e2] = eig(H_ff);  

e1 = diag(e1);  e2 = diag(e2); PP = blkdiag(1,P_e,P_f);

PP_full = blkdiag(eye(Nv),kron(P_e,eye(Nv)),kron(P_f,eye(Nv)));
%projectors in site only basis for single and excited states
alpha_n = mtimesx(mtimesx(PP,'C',V_n),PP); 
Cg  = squeeze(alpha_n(1,2:N+1,:)) ; %contract over dimension 2
Cf = alpha_n(2:N+1,N+2:end,:); %contract over dimension 3

 %diagonalise the blocks of the full Hamiltonian
[M_g,E0] = eig(H_vib); [M_e,E1] = eig(H_e);  [M_f,E2] = eig(H_f); 

E0 = diag(E0);  E1 = diag(E1);  E2 = diag(E2); MM = blkdiag(M_g,M_e,M_f);

alpha_n2 = mtimesx(mtimesx(MM,'C',V_n2),MM); 
Cg2  = alpha_n2(1:Nv,Nv+1:(N+1)*Nv,:) ; %contract over dimension 3
Cf2 = alpha_n2(Nv+1:(N+1)*Nv,(N+1)*Nv+1:end,:); %contract over dimension 3

%coefficients mapping the interaction operator into the new basis

lam_ex_full = zeros(1+N+N*(N-1)/2,1); lam_tot = zeros(N,1);
for k = 1:N
    tmp = diag(mtimesx(MM,'C',mtimesx(diag(fock_space_rep(:,k)),MM)));
   lam_ex_full = lam_ex_full + tmp.*lam_tot(k);
end

%% Calculate dipole moments averages in the exciton basis and such like
%1st order averages
av_set_fo = mtimesx(mtimesx(Cg,av_2),Cg,'C'); 
av_set_cd = mtimesx(mtimesx(Cg,av_3+av2_3),Cg,'C'); 

%also for the ex-vib basis transitions
Cg2_flat = reshape(Cg2,Nv^2*N,N); 
Av_set_fo = mtimesx(mtimesx(Cg2_flat,av_2),Cg2_flat,'C'); 
Av_set_cd = mtimesx(mtimesx(Cg2_flat,av_3+av2_3),Cg2_flat,'C'); 
% Capital letters denote the dipole averages in the full exciton
% vibrational basis

% How to proceed now actually matters what basis we use due to
% approximations used about the time evolution during t2

%3rd order averages
if strcmp(opts(1),'ex_basis') %compute in the exciton basis
tmp = mtimesx(mtimesx(Cg,av_4),Cg,'C'); %first two transitions are between
%g and one exciton transitions
tmp = permute(tmp,[3,4,1,2,5]); %permute so next two dimensions are in mtimesx

%the first one involves transitions between g-e only (GSB, SE)
av_set_GSB = mtimesx(mtimesx(Cg,squeeze(tmp)),Cg,'C');
av_set_GSB = permute(av_set_GSB,[3,4,1,2,5]); %permute back to original order

Cf_flat = reshape(Cf,N^2*(N-1)/2,N); %flatten out for use with mtimesx
%new order like Cf_flat(k,f,j) = Cf_flat(k+N*(f-1),j)

av_set_ESA = mtimesx(mtimesx(Cf_flat,squeeze(tmp)),Cf_flat,'C'); %ESA type
av_set_ESA = permute(av_set_ESA,[3,4,1,2,5]);

tmp = mtimesx(mtimesx(Cg,av_5+av2_5),Cg,'C');
tmp = permute(tmp,[3,4,1,2,5,6]); 

%approximate amplitude of k1, k2 by the vertical transition energy of the
%exciton states this map into k~w_res/c
%next I want to really scale SE and GSB differently to account for Stokes
%shift, this is a closer approximation Do this for each path when
%performing the loop
% E_GSB = E1+lam_ex_full(2:N+1);  
% E_SE = E1-lam_ex_full(2:N+1);  
av_set2_GSB  = mtimesx(mtimesx(Cg,tmp), Cg,'C');
av_set2_GSB = permute(av_set2_GSB,[3,4,1,2,5,6]);
av_set2_ESA  = mtimesx(mtimesx(Cf_flat,tmp), Cf_flat,'C');  
av_set2_ESA = permute(av_set2_ESA,[3,4,1,2,5,6]);


elseif strcmp(opts(1),'ex_vib_basis') %exciton vibrational basis

tmp2 = mtimesx(mtimesx(Cg2_flat,av_4),Cg2_flat,'C');
tmp2 = permute(tmp2,[3,4,1,2,5]); %permute so next two dimensions are in mtimesx

%the first one involves transitions between g-e only (GSB, SE)
av_set_GSB = mtimesx(mtimesx(Cg2_flat,squeeze(tmp2)),Cg2_flat,'C');
av_set_GSB = permute(av_set_GSB,[3,4,1,2,5]); %permute back to original order

%the next loop is more complicated because it involves transitions from e-f
%and then f'-e', need to reshape

Cf2_flat = reshape(Cf2,Nv^2*N^2*(N-1)/2,N); %flatten out

av_set_ESA = mtimesx(mtimesx(Cf2_flat,squeeze(tmp2)),Cf2_flat,'C'); %ESA type
av_set_ESA = permute(av_set_ESA,[3,4,1,2,5]);
%av_set_ESA (e1,e2,e3*(N-1)*N/2+f1,e4*(N-1)*N/2+f2) type calls
%next calculate the higher order moments from the other interactions

tmp2 = mtimesx(mtimesx(Cg2_flat,av_5+av2_5),Cg2_flat,'C');
tmp2 = permute(tmp2,[3,4,1,2,5]);

av_set2_GSB  = mtimesx(mtimesx(Cg2_flat,tmp2), Cg2_flat,'C');
av_set2_GSB = permute(av_set2_GSB,[3,4,1,2,5,6]);
av_set2_ESA  = mtimesx(mtimesx(Cf2_flat,tmp2), Cf2_flat,'C');  
av_set2_ESA = permute(av_set2_ESA,[3,4,1,2,5,6]);
elseif strcmp(opts(1),'site_basis') %work out in the site basis, 
    %can only be paired with Forster or pure dephasing
    
av_set_GSB = av_4;
Cf_flat = reshape(Cf,N^2*(N-1)/2,N); %flatten out for use with mtimesx
av_set_ESA = mtimesx(mtimesx(Cf2_flat,av_4),Cf2_flat,'C');
av_set_ESA = permute(av_set_ESA,[3,4,1,2,5]);
    
av_set2_GSB = av_5+av2_5;
av_set2_ESA  = mtimesx(mtimesx(Cf_flat,av_5+av2_5), Cf_flat,'C');      
av_set2_ESA = permute(av_set2_ESA,[3,4,1,2,5,6]);    
    
else
    error('opts(1) should contain either "ex_basis" "ex_vib_basis" or "site_basis"')
end
%{
%1st order averages
av_set_fo = mtimesx(mtimesx(M_e,'T',av_2),M_e); 
av_set_cd = mtimesx(mtimesx(M_e,'T',av_3+av2_3),M_e); 

%3rd order averages

tmp = mtimesx(mtimesx(M_e,'T',av_4),M_e); %first two transitions are between
%g and one exciton transitions

sz = size(tmp); 
if ndims(tmp) ==4
    sz = [sz,1];
end

tmp = permute(tmp,[3,4,1,2,5]); %permute so next two dimensions are in mtimesx
av_set_GSB = mtimesx(mtimesx(M_e,'T',tmp),M_e);
av_set_GSB = permute(av_set_GSB,[3,4,1,2,5]); %permute back to original order

av_set_ESA =  zeros(N,N,N,N,N*(N-1)/2,N*(N-1)/2,size(av_set_GSB,5));

for k3=1:N %much bette loop
    for k4=1:N
        tmp2 = zeros(N*(N-1)/2,N*(N-1)/2,sz(1),sz(2),sz(5));
        for flp = 1:N*(N-1)/2
            ss = find(fock2(flp,:));
            for flp2 = 1:N*(N-1)/2
                ss2 = find(fock2(flp2,:));
                
tmp2 = tmp2 + mtimesx(M_f(flp,:).'*M_f(flp2,:),(...
        M_e(ss(1),k3)*M_e(ss2(1),k4)*tmp(ss(2),ss2(2),:,:,:) +...
        M_e(ss(2),k3)*M_e(ss2(1),k4)*tmp(ss(1),ss2(2),:,:,:) +...
        M_e(ss(1),k3)*M_e(ss2(2),k4)*tmp(ss(2),ss2(1),:,:,:) +...
        M_e(ss(2),k3)*M_e(ss2(2),k4)*tmp(ss(1),ss2(1),:,:,:)));

            end
        end
        av_set_ESA(:,:,k3,k4,:,:,:) = permute(tmp2,[3,4,1,2,5]);        
    end
end

tmp = mtimesx(mtimesx(M_e,'T',av_5+av2_5),M_e);
tmp = permute(tmp,[3,4,1,2,5,6]); 

av_set2_GSB = mtimesx(mtimesx(M_e,'T',tmp),M_e);
av_set2_GSB = permute(av_set2_GSB,[3,4,1,2,5,6]);

sz = size(av_set2_GSB);
av_set2_ESA =  zeros(N,N,N,N,N*(N-1)/2,N*(N-1)/2,sz(5),sz(6));

for k3=1:N 
    for k4=1:N
        tmp2 = zeros(N*(N-1)/2,N*(N-1)/2,sz(1),sz(2),sz(5),sz(6));
        for flp = 1:N*(N-1)/2
            ss = find(fock2(flp,:));
            for flp2 = 1:N*(N-1)/2
                ss2 = find(fock2(flp2,:));
                
tmp2 = tmp2 + mtimesx(M_f(flp,:).'*M_f(flp2,:),(...
        M_e(ss(1),k3)*M_e(ss2(1),k4)*tmp(ss(2),ss2(2),:,:,:,:) +...
        M_e(ss(2),k3)*M_e(ss2(1),k4)*tmp(ss(1),ss2(2),:,:,:,:) +...
        M_e(ss(1),k3)*M_e(ss2(2),k4)*tmp(ss(2),ss2(1),:,:,:,:) +...
        M_e(ss(2),k3)*M_e(ss2(2),k4)*tmp(ss(1),ss2(1),:,:,:,:)));

            end
        end
        av_set2_ESA(:,:,k3,k4,:,:,:) = permute(tmp2,[3,4,1,2,5,6]);        
    end
end
%}
% Next Calculate 2D signals

%% Calculate operator propogating in time
H_tot = blkdiag(H_vib,H_e,H_f);
if strcmp(opts(1),'site_basis') 
    H_ex = H_tot;
L_full =  -1i*(kron(eye(length(H_tot)),H_tot)-kron(H_tot.',eye(length(H_tot))));
Lindblad_op = Lindblad_op_gen(B,BO_modes,[],[],1+N+N*(N-1)/2);    
elseif strcmp(opts(1),'ex_basis')
 H_ex =  PP_full'*H_tot*PP_full;                          
L_full =  -1i*(kron(eye(length(H_ex)),H_ex)-kron(H_ex.',eye(length(H_ex))));
Lindblad_op = Lindblad_op_gen(B,BO_modes,[],[],1+N+N*(N-1)/2,PP_full);   
elseif strcmp(opts(1),'ex_vib_basis')
  H_ex = diag([E0;E1;E2]);                       
L_full =  -1i*(kron(eye(length(H_ex)),H_ex)-kron(H_ex.',eye(length(H_ex))));
Lindblad_op = Lindblad_op_gen(B,BO_modes,[],[],1+N+N*(N-1)/2,MM);     
end


if strcmp(opts(2),'HEOM')
% calculate operator coupling excitations on each site to the exciton
%hamiltonian

V_coup_tot = zeros(Nv^2*(1+N+N*(N-1)/2),Nv^2*(1+N+N*(N-1)/2),N);
    
for j = 1:N %project each term into the appropriate basis
    if strcmp(opts(1),'ex_vib_basis')
        V_1 = M_e'*kron(diag(fock1(:,j)),eye(Nv))*M_e;
        V_2 = M_f'*kron(diag(fock2(:,j)),eye(Nv))*M_f;      
    elseif strcmp(opts(1),'ex_basis')
        V_1 = kron(P_e'*diag(fock1(:,j))*P_e,eye(Nv));
        V_2 = kron(P_f'*diag(fock2(:,j))*P_f,eye(Nv));
    else %must be site_basis
        V_1 = kron(diag(fock1(:,j)),eye(Nv));
        V_2 = kron(diag(fock2(:,j)),eye(Nv));   
    end
V_coup_tot(:,:,j) = blkdiag(0,V_1,V_2);  

end            

prop_op_full = H_prop_gen2(H_ex,H_ex,V_coup_tot,V_coup_tot,QQ...
                    ,const_factor,coup_com_save,coup_acom_save); 
prop_op_full = prop_op_full + kron(speye(HL),Lindblad_op);     

elseif strcmp(opts(2),'Redfield') 
if strcmp(opts(1),'site_basis')
    error('site_basis and Redfield cannot be used together')
elseif strcmp(opts(1),'ex_vib_basis')  

  [~,R_red_op] = redfield_calc_2(H_e,blkdiag(H_vib,H_f),...
            fock_space_rep,B,Drude_modes,[],BO_modes,[],[],true);   
    
elseif strcmp(opts(1),'ex_basis')      
sz1 =  1 + length(H_ee) + length(H_ff);
Nv =  length(H_vib);

R_red = redfield_calc_2(H_ee,blkdiag(0,H_ff),...
            fock_space_rep,B,Drude_modes,[],BO_modes,[],[],true);   
LL = sz1*Nv;
R_red_op = zeros(LL^2);

for a = 1:sz1
    for b = 1:sz1

    temp = squeeze(R_red(a,b,:,:));
    
        for av = 1:Nv
            for bv = 1:Nv
                temp2 = sparse(av,bv,1,Nv,Nv); %this vibrational element
                temp3 = kron(temp,temp2); %combine vib and 
                
    lp = ((b-1)*Nv + bv-1)*sz1*Nv + (a-1)*Nv + av; %this element is mapped to

    R_red_op(lp,:) = reshape(temp3,[LL^2,1]).';
            end
        end
    end
end                      
end
prop_op_full = L_full + Lindblad_op - R_red_op;
    
elseif strcmp(opts,'Pure_dephasing') 
        %includes a [V,[V,rho]] style term, still given the same variable
        %name
        
        R_red_op = sparse([],[],[],sz_full^2,sz_full^2);
        
      for j = 1:length(Drude_modes)
          
          QQ_j = QQ{j}; %QQ will be two numbers, one with a
          %pure dephasing part and another with a dissipative coefficient
          %site basis interaction operator to the exciton basis
             V = sparse(MM_full*kron(diag(fock_space_rep(:,j)),speye(Nv))*MM_full');
           Qcom =  kron(speye(size(V)),V)-kron(V.',speye(size(V)));
           Qacom = kron(speye(size(V)),V)+kron(V.',speye(size(V))) ; %anti commutator
           
            R_red_op = R_red_op + QQ_j(1).*Qcom*Qcom + QQ_j(2).*Qacom*Qcom;
          
      end      
      R_red_op = sparse(R_red_op);
      prop_op_full = L_full + Lindblad_op - R_red_op;
      
%elseif strcmp(opts,'Forster')       %maybe write some shit at some point
%elseif strcmp(opts,'Modified_Redfield')       %maybe write some shit at some point
end
%% Select slices of the operator for different blocks

%prop_op_full = L_full + Lindblad_op - R_red_op2;

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

%% Calculate the amplitude of the window functions
% Due to the frequency resolved detection we assume we have a delta
% function pulse in time for the probe but also full frequency resolution,
% I.e. E(t)E(t+t_3) = delta(t-t_sep) e^(i omega_3 t_3)   hence
% iint dt dt_3 E(t)E(t+t_3) [V_ge(t_3) V_eg] G(t) = V_ge_f(omega_3) V_eg
tic 
[V_ge_f,V_ef_f] = window_fun_HEOM_f3(N,HL,conj(permute(V_eg,[2,1,3])),...
                    conj(permute(V_fe,[2,1,3])),om3_rng,prop_op_full);
                clear prop_op_full
                
time_dep_of_output_op_calc_time = toc;
if t_verbose; time_dep_of_output_op_calc_time   %#ok<NOPRT>
end
%This is expressed in Liouville space                                

%also use this time dependent value to calculate linear spectra 
% dt3 =  t3_range(2)-t3_range(1); om3 = pi*(-1/dt3:2/t3_range(end):1/dt3);

sigma_om3 = zeros(size(V_ge_f,1),1); alpha_om3 = zeros(size(V_ef_f,1),1);
for j=1:size(V_ge_f,3)

        rho_sec = sparse(V_eg(:,:,j)*rho_0);
        rho_sec = reshape(rho_sec,numel(rho_sec),1);

    for j2 = 1:size(V_ge_f,3)
        sigma_om3 = sigma_om3 + av_set_fo(j,j2)*V_ge_f(:,1:N*Nv^2,j2)*rho_sec;
        alpha_om3 = alpha_om3 + av_set_cd(j,j2)*V_ge_f(:,1:N*Nv^2,j2)*rho_sec;
    end
end
alpha_om3 = alpha_om3.*reshape(om3_rng,size(alpha_om3));
lin_spec={sigma_om3,alpha_om3};
%% first calculate the initial coherences set up

%generate the N initial density matricies from the first pathways

if strcmp(opts(1),'ex_basis') || strcmp(opts(1),'site_basis') 
    %each vibrational state is included in the same transitions
    %operators are the same as diagonalisation in implicit in dipole
    %averaging
    V_eg = V_n2(j*Nv+1:(N+1)*Nv,1:Nv,:); 
    V_fe = V_n2((N+1)*Nv:end,j*Nv+1:(N+1)*Nv,:); 
elseif strcmp(opts(1),'ex_vib_basis') 
    %transitions go from exvib - exvib
    % Nv -> N*Nv transitions and N*Nv - >N*(N-1)/2*Nv transitions
    V_eg = zeros(N*Nv,Nv,Nv*N*Nv); V_fe = zeros(N*(N-1)*Nv/2,N*Nv,N*Nv*N*(N-1)*Nv/2);
%this is actually a pretty inefficient way to store these....
cntge = 0; cntef = 0;
    for lp1 = 1:Nv
        for lp2 = 1:N*Nv
            cntge  = cntge +1;
            V_eg(lp2,lp1,cntge) = 1;
        end
    end
    for lp1 = 1:N*Nv
        for lp2 = 1:N*(N-1)/2*Nv
            cntef  = cntef +1;
            V_fe(lp2,lp1,cntef) = 1;
        end
    end       
end
%{
for j=1:N
    tmp = V_n2(j*Nv+1:(N+1)*Nv,1:Nv); 
    tmp2 = V_n2((N+1)*Nv:end,j*Nv+1:(N+1)*Nv); 
    if strcmp(opts(1),'ex_vib_basis')      
        V_eg(:,:,j) = M_e'*tmp*M_g;
        V_fe(:,:,j) = M_e'*tmp*M_g;
    elseif strcmp(opts(1),'ex_basis')
        V_eg(:,:,j) = kron(P_e',eye(Nv))*tmp*kron(P_g,eye(Nv));
    else %must be site_basis
        V_eg(:,:,j) =tmp;
    end
end
%}
 tic
 
 rho_0 = expm(-B*Hvib); %thermal initial state

 if ~isempty(om_1_rng)  %prop in frequency space
[rho_eg_f,rho_ge_f] = rho_eg_f_HEOM3(N,Nv,V_eg,rho_0,om1_rng...
                                        ,prop_op_eg,prop_op_ge,HL);
 else %need to make this the correct empty size
    rho_eg_f = zeros(N*Nv^2*HL,0,N); rho_ge_f = rho_eg_f;
 end
 if ~isempty(t1_range) %prop in time space (relevant for short pulses)
 [rho_eg_t,rho_ge_t] = rho_eg_t_HEOM3(N,Nv,V_eg,rho_0,...
                                t1_range,prop_op_eg,prop_op_ge,HL); 
 else %this will be empty
    rho_eg_t = zeros(N*Nv^2*HL,0,N); rho_ge_t = rho_eg_t;                         
 end
 %combine together, frequency then time space.  Generally this is expected
 %to be used as "either" prop in time or freq OR just the zero value for
 %time used to calculate the impulsive response.
 rho_eg = [rho_eg_f,rho_eg_t];  rho_ge = [rho_ge_f,rho_ge_t];
 
end
init_coherence_time_to_calc = toc;
if t_verbose
    init_coherence_time_to_calc %#ok<NOPRT>
end


%% Now calculate each of the required neq second order matricies

%this method works by calculating slices of the density matrix over t1 (or
%omega1 in frequency space) 

sz_params = [N,Nv,HL]; %number of sites, vibrations length of HEOM

[GSB_s,SE_s,ESA_s,GSB_c,SE_c,ESA_c,rho_gg_t,rho_ee_t] = PP_cont...
         (sz_params,rho_eg,rho_ge,V_eg,V_ef,prop_op_gg,prop_op_ee,t2_range...
         ,E1,E2,av_set_GSB,av_set_ESA,av_set2_GSB,av_set2_ESA,HL,[],t_verbose) ;
   

%% Add together to create average

if ~exist('GSB_sav','var') %preallocate as it is the first time
    
    GSB_sav = GSB_s;  SE_sav = SE_s;   ESA_sav = ESA_s;   
    GSB_cav = GSB_c;  SE_cav = SE_c;   ESA_cav = ESA_c;   
else %add extra terms
    GSB_sav = GSB_sav+GSB_s;  SE_sav = SE_sav+SE_s;   ESA_sav = ESA_sav+ESA_s;   
    GSB_cav = GSB_cav+GSB_c;  SE_cav = SE_cav+SE_c;   ESA_cav = ESA_cav+ESA_c;   
end
end