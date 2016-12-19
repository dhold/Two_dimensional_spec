function [Pol_t,rho_g,rho_e,rho_f,rho_t] = non_perturbative_spec...
                (Temp,beam_param,sys_params,t_range_fs,phases,tau1,tau2,rot_mat...
                 ,max_tier,Kappa, useRWA,om_s,sd_mat,tol,clever_prop)
%calculates the polarization (and third order density matrix) at a
%particular point, determined by the phases which are input. These
%are related by  phases = [phi_1 + k1.r ,phi_2 + k2.r,..] etc 
% tau1 and tau2 are delays between beam 1 and 2 and 2 and 3  in fs
% Kappa is the # of matsubara frequencies max_tier is max Heirarchy tier
%useRWA whether or not to use the rotating wave approx, om_s is a scaling
%frequency, sd_mat is a shift to energies and tol is the solver tolerance
if ~exist('clever_prop','var')
clever_prop = false; 
end
%Temp   %temperature in Kelvin
[t_sc, B,units,unitnames]= inv_cm_unit_sys(Temp); 
E_scale = units{2}^2*units{3}/units{4}^2; %energy scaling

%set of 1 by 3 vectors which control the beam parameters
FWHM = beam_param{1}; %full width at half max of pulse intensity profile in fs
Peak_intensity = beam_param{2}; %stupid units of MW cm^{-2} nm^{-1}
om_carrier = beam_param{3}; om = om_carrier;%carrier frequencies
pol = beam_param{4}; %polarization of the beams

Lambda= 1e7./om_carrier; %wavelength in nm
KK = 0.441; %some shape dependent constant
delta_lambda = (Lambda).^2/(FWHM*1e-6*units{1}/2/pi)*KK; %bandwidth (nm)

sig = FWHM/(4*sqrt(log(2)))*t_sc*1e-3;  %Gaussian param for E in inverse cm
%|peak E field| e * L_0 = 2 sqrt(alpha / n)*sqrt(h * (I* delta lam *L_0^2))
% L_0 = 1 A  = 1e-8 cm
E0 = 2*sqrt(7.2973525698e-3/1.33).*... %units inverse cm 
    sqrt(units{5}.*(1e6*Peak_intensity.*1e-16).*(delta_lambda))/E_scale; 

t_range = t_range_fs * t_sc *1e-3; 

%%

%note if E.mu ~ H_site coupling is strong

%frequencies will ultimately be scaled by om_s where required, a better
%choice than the mean carrier frequency can be found

%electric field envelopes

%scale into inverse cm unit system
tau1 = tau1 * t_sc *1e-3; 
tau2 = tau2 * t_sc *1e-3; 

%Assumed envelope real, if you change this change the operators to account
%for it
E1env = @(t) E0(1)*exp(-t.^2/sig(1)^2/2);
E2env = @(t) E0(2)*exp(-(t-tau1).^2/sig(2)^2/2);
E3env = @(t) E0(3)*exp(-(t-tau2-tau1).^2/sig(3)^2/2);

if ischar(sys_params)
load(sys_params,'H_site','dipole_data')
load(sys_params,'Drude_modes','BO_modes','pdm_shift')
else
H_site = sys_params{1};     dipole_data = sys_params{2};  
Drude_modes = sys_params{3};  BO_modes = sys_params{4};  
pdm_shift = sys_params{5}; 
end

N = length(H_site); %number of sites
if isempty(rot_mat)
   rot_mat = eye(3); %identity rotation 
end
R = rot_mat*dipole_data{4}.'; % Angstroms
mu = rot_mat*dipole_data{3}.'; % e * Angstroms  
mdip = cross(mu,R);  

%% HEOM stuff

clear cc vv cc2 cc_com cc_acom

g_cnst = zeros(N,1); lam_tot = zeros(N,1);
    for j = 1:N
        if ~isempty(BO_modes{j})
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
   (lambda,gam,om_0,Temp,Kappa,gam_dru,lam_dru); 

     cc_com{j}= [cc2R,cc1]; vv{j} = [vv2,vv1];%#ok<*AGROW> %cc(j,:)=cc_tmp; 
     cc_acom{j}= [cc2I,cc1*0]; QQ_topass(j,:) = [QQ,0]+QQ_extra;
     
    end
    
[coup_com_save,coup_acom_save,const_factor,nn]=...
    multi_site_drude_and_brownian_HOM_4(cc_com,cc_acom,vv,inf,max_tier);

%%  Full operators

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

fock1 = fock_space_rep(2:N+1,:); fock2 = fock_space_rep(N+2:end,:);
% Generate the operator components of the interaction dipole moment 
%associated with each transition
NN = 1+N+N*(N-1)/2;
V_n = zeros(NN,NN,N); V_coup_tot= V_n;
Vxyz =zeros(NN,NN,3); %full tensor operator

for lp =1:N
   
    V = zeros(NN); 
    V(lp+1,1) = 1;  %mixes to ground
    
    lg = [false(N+1,1);fock_space_rep(N+2:end,lp)==1]; %elements with excitation at lp
    lg2 = fock_space_rep; lg2(:,lp) = 0; %can't mix to itself 
    V(lg,2:N+1) =lg2(lg,:);
    
    Vxyz = Vxyz + cat(3,V.*mu(1,lp),V.*mu(2,lp),V.*mu(3,lp)); 
    
    V = V+V'; V_n(:,:,lp) = V;
    
    V_coup_tot(:,:,lp) = blkdiag(0,diag(fock1(:,lp)),diag(fock2(:,lp)));

end

H_g = 0; 
H_e = H_site +sd_mat;% add static disorder
H_f = zeros(N*(N-1)/2); %double excited Hamiltonian
occ_f = logical(fock_space_rep(sum(double(fock_space_rep),2)==2,:));

for k = 1:N*(N-1)/2 %enumerate the coupling Hamiltonian
    
        lg1 = occ_f(k,:); s1 = find(lg1);
       H_f(k,k)= H_e(s1(1),s1(1))+H_e(s1(2),s1(2)); 
       if ~isempty(pdm_shift) %shift to double excited states
           H_f(k,k) = H_f(k,k) + pdm_shift{s1(1),s1(2)};
       end
       
       for kk = k+1:N*(N-1)/2
           lg2 = occ_f(kk,:); lg3 = lg1 & lg2; %repeated element
            if any(lg3)
                s2 = find(lg2 &~lg3); s3 = find(lg1 &~lg3);
                H_f(k,kk) = H_e(s2,s3);
                H_f(kk,k) = H_f(k,kk)';
            end               
       end  
end
 H_tot = blkdiag(H_g,H_e,H_f); 

[M_e,E1] = eig(H_e);  [M_f,E2] = eig(H_f);  %diagonalise the blocks

%% generate operator for the time independent part

prop_op_full = H_prop_gen2(H_tot,H_tot,V_coup_tot,V_coup_tot,QQ_topass,...
                    const_factor,coup_com_save,coup_acom_save);  
       
HL = length(const_factor);  


tmp1 = false(NN); tmp1a = tmp1; 
tmp2 = tmp1; tmp2a = tmp2;  tmp3 = tmp1; tmp3a = tmp1;

tmp1(1,1) = 1;  tmp1a(2:N+1,1) = true; tmp1b = tmp1a';
tmp2(2:N+1,2:N+1) = true; tmp2a(N+2:end,2:N+1) = true; tmp2b =tmp2a';
tmp3(N+2:end,N+2:end) = true; tmp3a(N+2:end,1) = true; tmp3b = tmp3a';

sel_eg = tmp1a; sel_ge = tmp1b; sel_fe = tmp2a; sel_ef = tmp2b; %used at the end

tmp1 = padded(tmp1,HL); tmp1a = padded(tmp1a,HL); tmp1b = padded(tmp1b,HL); 
tmp2 = padded(tmp2,HL); tmp2a = padded(tmp2a,HL); tmp2b = padded(tmp2b,HL);
tmp3 = padded(tmp3,HL); tmp3a = padded(tmp3a,HL); tmp3b = padded(tmp3b,HL);
%individual sections, when no electric field is present the are all
%independent and propogation can be done separately 

%choose what basis to use depending on the steady state of first 
%exciton manifold 
try 
    [V1,~] = eigs(prop_op_full(tmp2,tmp2),1,'sm');
VV1 = reshape(V1(1:N^2),N,N); VV1 = VV1 / trace(VV1); %complex prefactor
ipr =  inverse_participation_ratio(VV1);

if ipr < 1/N + (1-1/N)*2/5
    use_exciton_basis = false; %mostly localized
else
    use_exciton_basis = true; %mostly delocalized
    %obvious stuff inbetween exists, particularly when some sites have very
    %low iprs and some high but hey
   
        
for lp =1:N
   
    V_1 = M_e'*diag(fock1(:,lp))*M_e;
    V_2 = M_f'*diag(fock2(:,lp))*M_f;
    
    V_coup_tot(:,:,lp) = blkdiag(0,V_1,V_2);

end    
 H_ex = blkdiag(0,E1,E2); 
    prop_op_full = H_prop_gen2(H_ex,H_ex,V_coup_tot,V_coup_tot,QQ_topass,...
                    const_factor,coup_com_save,coup_acom_save);  
    
    
end
 catch ME
    use_exciton_basis = false;
end

if clever_prop
sep_prop_free({tmp1,tmp1a, tmp1b, tmp2,tmp2a, tmp2b,tmp3,tmp3a,tmp3b},...
                prop_op_full,om_s); %pass params
end

%%

%interaction operators

V1 = Vxyz(:,:,1)*pol(1,1)+Vxyz(:,:,2)*pol(1,2)+Vxyz(:,:,3)*pol(1,3); 
V2 = Vxyz(:,:,1)*pol(2,1)+Vxyz(:,:,2)*pol(2,2)+Vxyz(:,:,3)*pol(2,3);
V3 = Vxyz(:,:,1)*pol(3,1)+Vxyz(:,:,2)*pol(3,2)+Vxyz(:,:,3)*pol(3,3);

if use_exciton_basis
    %project to ze new basis
    PP = blkdiag(1,M_e,M_f);
    V1 = PP'*V1*PP;  V2 = PP'*V2*PP;  V3 = PP'*V3*PP;
end

padmat = speye(HL);

if ~useRWA %don't use RWA, slower

VV1 = -1i*sparse(kron((V1+V1').',speye(length(V1)))...
        -kron(speye(length(V1)),(V1+V1')) ); 
VV2 = -1i*sparse(kron((V2+V2').',speye(length(V2)))...
        -kron(speye(length(V2)),(V2+V2')) );     
VV3 = -1i*sparse(kron((V3+V3').',speye(length(V3)))...
        -kron(speye(length(V3)),(V3+V3')) );         
    
VV1 = kron(padmat, VV1); VV2 = kron(padmat, VV2);  VV3 = kron(padmat, VV3);      

%Assumed that envelope fun is real or additional sin terms

LI = @(t) (E1env(t).*cos(phases(1) - om(1)*t))*VV1  +...
        (E2env(t).*cos(phases(2) - om(2)*t))*VV2 +...
        (E3env(t).*cos(phases(3) - om(3)*t))*VV3; %actual non RWA time-dep lio

else %use RWA, should be fine if the fields are not super strong   
    
VV1f =  1i*sparse(kron(speye(length(V1)),V1)-kron(V1.',speye(length(V1))));
VV1b =  -1i*sparse(kron((V1').',speye(length(V1)))-kron(speye(length(V1)),V1'));
VV2f =  1i*sparse(kron(speye(length(V2)),V2)-kron(V2.',speye(length(V2))));
VV2b =  -1i*sparse(kron((V2').',speye(length(V2)))-kron(speye(length(V2)),V2'));
VV3f =  1i*sparse(kron(speye(length(V3)),V3)-kron(V3.',speye(length(V3))));
VV3b =  -1i*sparse(kron((V3').',speye(length(V3)))-kron(speye(length(V3)),V3'));


VV1f = kron(padmat, VV1f); VV1b = kron(padmat, VV1b);        
VV2f = kron(padmat, VV2f); VV2b = kron(padmat, VV2b);    
VV3f = kron(padmat, VV3f); VV3b = kron(padmat, VV3b); 
if 1==0
VV1f(17:end) = 0; VV1b(17:end) = 0;     
VV2f(17:end) = 0; VV2b(17:end) = 0;    
VV3f(17:end) = 0; VV3b(17:end) = 0;     
    
end
%Assumed envelope real

LI = @(t)   VV1f.*(E1env(t).*exp(+1i*(phases(1) - om(1)*t))/2)  +...
            VV1b.*(E1env(t).*exp(-1i*(phases(1) - om(1)*t))/2)  +...
            VV2f.*(E2env(t).*exp(+1i*(phases(2) - om(2)*t))/2)  +...
            VV2b.*(E2env(t).*exp(-1i*(phases(2) - om(2)*t))/2)  +...
            VV3f.*(E3env(t).*exp(+1i*(phases(3) - om(3)*t))/2)  +...
            VV3b.*(E3env(t).*exp(-1i*(phases(3) - om(3)*t))/2);       
%{
LI = @(t)   VV1f.*(E1env(t).*exp(+1i*(phases(1) - om(1)*t))/2)  +...
            VV1b.*(conj(E1env(t)).*exp(-1i*(phases(1) - om(1)*t))/2)  +...
            VV2f.*(E2env(t).*exp(+1i*(phases(2) - om(2)*t))/2)  +...
            VV2b.*(conj(E2env(t)).*exp(-1i*(phases(2) - om(2)*t))/2)  +...
            VV3f.*(E3env(t).*exp(+1i*(phases(3) - om(3)*t))/2)  +...
            VV3b.*(conj(E3env(t)).*exp(-1i*(phases(3) - om(3)*t))/2);           
%}         
end    
%calculate scaling vector M1 

Z0 = blkdiag(ones(length(H_g)),ones(length(H_e)),ones(length(H_f))); %non scaled

Z1= 0*Z0; Z3= Z1;  ZE = Z1; %scaled by different amounts
Z0 = reshape(Z0,numel(Z0),[]); Z0 = repmat(Z0,HL,1);

Z1(2:N+1,1) = 1; Z1(N+2:end,2:N+1) = 1; Z2 = Z1';  %scaled by w_s
Z3(N+2:end,1) = 1; Z4 = Z3'; %scaled by 2*w_s
Z1 = reshape(Z1,numel(Z1),[]); Z2 = reshape(Z2,numel(Z2),[]);
Z3 = reshape(Z3,numel(Z3),[]); Z4 = reshape(Z4,numel(Z4),[]);
Z1 = repmat(Z1,HL,1); Z2 = repmat(Z2,HL,1);
Z3 = repmat(Z3,HL,1); Z4 = repmat(Z4,HL,1);

% Scale the ODE so it is numerically easier to solve
%{
% M1 = @(t) Z0 + Z1*exp(+1i*om_s*t) + Z2*exp(-1i*om_s*t)...
%              + Z3*exp(+2i*om_s*t) + Z4*exp(-2i*om_s*t);
% M1inv = @(t) Z0 + Z1*exp(-1i*om_s*t) + Z2*exp(+1i*om_s*t)...
%              + Z3*exp(-2i*om_s*t) + Z4*exp(+2i*om_s*t);         
%                  %and it's derivative
% M2 = @(t) Z1*(1i*om_s*exp(+1i*om_s*t)) - Z2*(1i*om_s*exp(-1i*om_s*t))...
%         + Z3*(2i*om_s*exp(+2i*om_s*t)) - Z4*(2i*om_s*exp(-2i*om_s*t)); 
%    
%}
meth_of_scale = 1; %note this requires changing the code
 TAU = 2*pi; 
if meth_of_scale==1
   %slightly different way to calculate this as exp is stupid
M1 = @(t) Z0 + Z1*exp(+1i*mod(om_s*t,TAU)) + Z2*exp(-1i*mod(om_s*t,TAU))...
             + Z3*exp(+2i*mod(om_s*t,TAU)) + Z4*exp(-2i*mod(om_s*t,TAU));     
                 %and it's derivative
M2 = @(t) 1i*om_s*(Z1*(exp(+1i*mod(om_s*t,TAU))) - Z2*(exp(-1i*mod(om_s*t,TAU)))...
        + Z3*(2*exp(+1i*mod(2*om_s*t,TAU))) - Z4*(2*exp(-1i*mod(2*om_s*t,TAU))));     
    
%in terms of time difference so they don't scale stupidly at late times
%(not actually a problem it seems...)
%{
MM1 = @(t,Mold,told) Mold.*(Z0 + Z1*exp(+1i*om_s*(t-told)) + Z2*exp(-1i*om_s*(t-told))...
             + Z3*exp(+2i*om_s*(t-told)) + Z4*exp(-2i*om_s*(t-told)));
MM2 = @(t,Mold,told) 1i*om_s*Mold.*(Z1*exp(+1i*om_s*(t-told)) - Z2*exp(-1i*om_s*(t-told))...
        + Z3*2*exp(+2i*om_s*(t-told)) - Z4*2*exp(-2i*om_s*(t-told))); 
%}

% Can also scale coherences between the sites, this should speed stuff up
% unless the vibrational coupling is absurd, but will take longer to do the
% scaling stuf
elseif meth_of_scale == 2

warning('currently this isnt even written')

else  %test without scaling, should give the same populations but be slower
M1 = []; %    M1 = @(t) ones(size(Z0));
M2 = []; %    M2 = @(t) zeros(size(Z0));
end

scaled_ODE_fn({M1,M2,LI,prop_op_full}); %pass scaling matricies over

rho0 = zeros(length(prop_op_full),1); rho0(1) = 1; %simple initial condition
%simple method
%tic
if ~clever_prop %normal propogation the whole way
rho_t = OD_wrapper(t_range,@scaled_ODE_fn,rho0,[],'ode45',tol);%xnumel(H_sc));
rho_t = rho_t.'; 
%take range up to 7 standard deviation before first and 5 after last pulse
%e.g. t_range = linspace(-7*sig1,tau+5*sig2,2);

% Calculate polarization, requires only tier zero component of HEOM

Pol_t = zeros(3,length(t_range)); %note that the rescaling must be undone

if ~isempty(M1)
for k = 1:length(t_range); %undo the rescaling on each element
            tmp = M1(t_range(k)); 
            rho_t(:,k) = rho_t(:,k)./tmp;  
end
end

for j=1:3
    
    tmp_op = reshape(Vxyz(:,:,j)+Vxyz(:,:,j)',[NN^2,1])';
   %tmp_op = reshape(Vxyz(:,:,j),[NN^2,1])';
    Pol_t(j,:) = tmp_op*rho_t(1:NN^2,:);
    %This is just the total polarization at a point (which is effectively
    %determined by the relative phases between 
    
end
%toc
if nargout >1 
    rho_g = rho_t(tmp1,:); rho_e = rho_t(tmp2,:); rho_f = rho_t(tmp3,:);
end
    
else
%more complex method using the fact that when the E field is weak we can
%propogate the different components seperately and only prop the coherences
%at the end of the time range

np = floor(length(t_range)/2)+1; %number of points
if E0(1) ~=0
t1_range = [min(t_range(1),-7*sig(1)),min(6*sig(1),t_range(end))];

rho_t1 = OD_wrapper({t1_range,np},@scaled_ODE_fn,rho0,[],'ode45',tol);%xnumel(H_sc));
t1 = rho_t1{1}; rho_t1 = rho_t1{2};  %t1e = t1(end) 

%now propogate each component seperately between pulses, assuming this gap
%exists
else
    t1 = min(t_range(1),-7*sig(1));
    rho_t1 = rho0.';
end
t1_free = [t1(end),min(t_range(end),min(tau1-6*sig(2),tau1+tau2 -6*sig(3)))];
    if t1_free(2) > t1_free(1) %some free propogation time

[rho_t1_free,t1f] = sep_prop_free(t1_free(1),t1_free(end),np,rho_t1(end,:).',tol);
t1f = reshape(t1f,[length(t1f),1]);
rho_t1 = [rho_t1(1:end-1,:);rho_t1_free]; t1 = [t1(1:end-1,:);t1f]; %t1fe = t1f(end)   
    end
    
if E0(2) ~=0
 t2_range = [t1(end),min(tau1+6*sig(2),t_range(end))];
    if t2_range(2) > t2_range(1) %prop over the second pulse time

rho_t2 = OD_wrapper({t2_range,np},@scaled_ODE_fn,rho_t1(end,:).',[],'ode45',tol);
t2 = rho_t2{1}; rho_t2 = rho_t2{2}; 
rho_t1 = [rho_t1(1:end-1,:);rho_t2];   t1 =[t1(1:end-1);t2]; %t2e = t2(end)   
    end
end

 t2_free = [t1(end),min(t_range(end),tau1+tau2 -6*sig(3))];
     if t2_free(2) > t2_free(1) %some free propogation time

[rho_t2_free,t2f] = sep_prop_free(t2_free(1),t2_free(end),np,rho_t1(end,:).',tol);
 t2f = reshape(t2f,[length(t2f),1]);
rho_t1 = [rho_t1(1:end-1,:);rho_t2_free]; t1 = [t1(1:end-1,:);t2f];   %t2fe =  t2f(end)    
     end
     
if E0(3) ~=0
 t3_range = [t1(end),min(t_range(end),tau1+tau2 +6*sig(3))];
     if t3_range(2) > t3_range(1) %prop over the third pulse time
rho_t3 = OD_wrapper({t3_range,np},@scaled_ODE_fn,rho_t1(end,:).',[],'ode45',tol);
t3 = rho_t3{1}; rho_t3 = rho_t3{2}; 
rho_t1 = [rho_t1(1:end-1,:);rho_t3];   t1 =[t1(1:end-1);t3];  %t3e =  t3(end) 
     end
end
     %finally prop over the rest of the time
     t3_free = [t1(end),t_range(end)];
     if t3_free(2) > t3_free(1) %some free propogation time
[rho_t3_free,t3f] = sep_prop_free(t3_free(1),t3_free(end),np,rho_t1(end,:).',tol);
t3f = reshape(t3f,[length(t3f),1]);
rho_t1 = [rho_t1(1:end-1,:);rho_t3_free];   t1 =[t1(1:end-1,:);t3f]; %t3fe =  t3f(end)  
     end       
    
clear rho_t2 rho_t2_free rho_t3 rho_t3_free     
rho_t = interp1(t1,rho_t1,t_range); %interp1(t1,rho_t1,t_range,'cubic');
clear rho_t1
rho_t = rho_t.'; %put time dimension second


%select the elements corresponding to coherences to calculate the output
%polarization

%sel_eg = tmp1a; sel_ge = tmp1b; sel_fe = tmp2a; sel_ef = tmp2b; 

Vlio = squeeze(reshape(Vxyz+permute(conj(Vxyz),[2,1,3]),[NN^2,3]))';
if 1==0
[Pol_t1,Pol_t2] = deal(zeros(3,length(t_range))); 
%don't rescale the coherences and given them seperately 
for j=1:3
    
    Pol_t1(j,:) = (Vlio(j,sel_eg)*rho_t(sel_ge,:) +...
                Vlio(j,sel_fe)*rho_t(sel_ef,:)); %fwd
    Pol_t2(j,:) = (Vlio(j,sel_ge)*rho_t(sel_eg,:) + ...
                Vlio(j,sel_ef)*rho_t(sel_fe,:) ); %backward term
    %This is just the total polarization at a point in space, I may also
    %need to consider triple excited states if required!
    
end
Pol_t = {Pol_t2,Pol_t1};
else
Pol_t = zeros(3,length(t_range)); 
%note that the rescaling must be undone

rs_fct = exp(1i*t_range*om_s); %scaling factor of the coherences

for j=1:3
    
    Pol_t(j,:) = rs_fct.*(Vlio(j,sel_eg)*rho_t(sel_ge,:) +...
                Vlio(j,sel_fe)*rho_t(sel_ef,:))+...
                conj(rs_fct).*(Vlio(j,sel_ge)*rho_t(sel_eg,:) + ...
                Vlio(j,sel_ef)*rho_t(sel_fe,:) );
    %This is just the total polarization at a point in space, I may also
    %need to consider triple excited states if required!
    
end
    
end


%toc
if nargout >1 
    rho_g = rho_t(tmp1,:); rho_e = rho_t(tmp2,:); rho_f = rho_t(tmp3,:);
end


end
end
function padout = padded(AA,BB)

    AA = reshape(AA,[numel(AA),1]);
    padout = repmat(AA,[BB,1]); 

end