function [fock_space_rep,H_el,H_vib,H_e,H_f,site_shift,...
    V_n,V2_n,av_2,av_3,av2_3,av_4,av_5,av2_5,Drude_modes,BO_modes] = ...
    TwoD_with_Redfield_prerun(system_parameters_file,...
    num_realisations,param_set,numvib) ; %#ok<*STOUT>



%% load system parameters from file
load(system_parameters_file,'H_site','sd_mat','BO_modes','Drude_modes',...
    'mu','R','mdip','mu_sd','pdm_shift'); %loads parameters
N = length(H_site); %number of sites
if exist('numvib','var') %set number of vibrations in each modes, modes with 
    %"0" in 4th value contribute only to incoherent processes
    for j =1:length(numvib)
       BO_modes{j}(:,4) = numvib{j};
    end
end

if all(mu_sd==0) %#ok<*NODEF>
    mu_sd = []; end

sd_noise = randn(num_realisations,length(sd_mat)); 
site_shift = sd_noise.*repmat(sd_mat,[num_realisations,1]); %sd_mat is diagonal uncorrelated noise
if ~isempty(mu_sd) %static disorder in mu,  assumed uncorrelated 
    %to site static disorder, stupid, but what can you do?
    warning('currently static disorder in mu is not supported')
    sd_noise2 = randn(num_realisations,length(sd_mat)); 
    mu_shift = sd_noise2.*repmat(mu_sd,[num_realisations,1]);
end

%generate shifts to energies of each double excited states if not saved
if ~exist('pdm_shift','var') 
   load(system_parameters_file,'pdm'); %might have pdms listed instead
   if exist('pdm','var') %perm dipole moments known but not shifts to states
        load(system_parameters_file,'R');
        for j =1:N
            for kk = 1:N
pdm_shift{kk,j} = dot(pdm(kk,:),pdm(j,:))-3*dot(pdm(kk,:),R(j,:))...
              *dot(pdm(j,:),R(kk,:))/norm(R(j,:)-R(kk,:))^2;
pdm_shift{kk,j} = pdm_shift/norm(R(j,:)-R(kk,:))^3  ;    %#ok<*AGROW>
            end
        end       %loop over states
   else
       pdm_shift=[]; %no perm dipole moment
   end
end

%% Generate vibrational Hamiltonian and coupling terms to the electronic degrees of freedom

[H_el,H_vib,H_e,H_f,~,~,fock_space_rep] = ...
    generate_ex_vib_ham3(H_site,BO_modes,[],[],[],[],[],pdm_shift) ;
   
%to include state disorder add
%kron(eye(size(H_vib),fock_space_rep * site_shift))

%reshape(P' * C * P, NN_1*NN_2,1) = kron(P.',P')*reshape(C, NN_1*NN_2,1)
%reshape(P * C * P', NN_1*NN_2,1) = kron((P').',P)*reshape(C, NN_1*NN_2,1)
% Will need to project this into the interaction basis

%Lindblad_op = Lindblad_op_gen(B,BO_modes,[],[],1+N+N*(N-1)/2);

%% Generate the operator components of the interaction dipole moment 
%associated with each transition
V_n = zeros(1+N+N*(N-1)/2,1+N+N*(N-1)/2,N); sz2 = length(H_vib);
V2_n = zeros(sz2*(1+N+N*(N-1)/2),sz2*(1+N+N*(N-1)/2),N);

for lp =1:N
   
    V = zeros(1+N+N*(N-1)/2);  
    V(1,lp+1) = 1;  %mixes to ground
    
    lg = [false(N+1,1);fock_space_rep(N+2:end,lp)==1]; %elements with excitation at lp
    lg2 = fock_space_rep; lg2(:,lp) = 0; %[~,lg2] = find(lg2(lg,:)); 
    V(lg,2:N+1) =lg2(lg,:);
    V = V+V'; V_n(:,:,lp) = V;
 
    V2 = kron(V,eye(size(H_vib)));
    V2_n(:,:,lp) = V2;
    
end
%% Precompute all the possible dipole expectation values in site basis
% of the form <mu_4 dot e_out* ... mu_1 dot e_1  exp(i*k dot r_4..)>_{isotropic}
% these are computed as taylor expansions exp(1ix) ~ 1 + ix in terms of 4th
% and 5th order tensor averages

if ~iscell(param_set)
[~,~,~,av_4,av_5,av2_5]=ori_precomp_site_basis(param_set,mu,R,imag(mdip));
else %compute different order with different param_set
[~,~,~,av_4,av_5,av2_5]=ori_precomp_site_basis2(param_set{1},param_set{2},mu,R,imag(mdip));
end
%also compute some first order averages with CP light
 tmp = [0,0,1;[1,+1i,0]/sqrt(2);[1,-1i,0]/sqrt(2)];
 
[av_2,av_3,av2_3]=ori_precomp_site_basis(tmp,mu,R,imag(mdip));
 

