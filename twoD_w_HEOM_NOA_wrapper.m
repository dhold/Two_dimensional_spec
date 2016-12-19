% This code calculates 2D spectroscopy using the HEOM in the exciton basis
% for an excitonic system with a spectral density of underdamped and
% overdamped brownian oscillator modes on each side.
% NOA stands for no orientation averaging, this version works for a
% particular orientation and a set of polarizations.  If you want
% orientation averages passing no polarization will calculate all the terms
% required for the most general set of averages.  If "bb" is given as
% a logical true value it will also calculate all the relevant terms for
% the magnetic field.  

%list of all required polarizations to get complete orientation averages
%for ANY configuration
clear pol
bb= true;
 x =[1,0,0]; y=[0,1,0]; z =[0,0,1];
required_pol = cat(3,[x;x;y;y],[y;y;x;x],[x;x;z;z],[z;z;x;x],[y;y;z;z],[z;z;y;y]);
required_pol = cat(3,[x;x;x;x],[y;y;y;y],[z;z;z;z],required_pol,...
                required_pol([1,3,2,4],:,:),required_pol([1,4,3,2],:,:));
 
%required_pol = {[x;x;y;y],[y;y;x;x],[x;x;z;z],[z;z;x;x],[y;y;z;z],[z;z;y;y],...
%           [x;y;x;y],[y;x;y;x],[x;z;x;z],[z;x;z;x],[y;z;y;z],[z;y;z;y],...     }
%           [x;y;y;x],[y;x;x;y],[x;z;z;x],[z;x;x;z],[y;z;z;y],[z;y;y;z],...
%            [x;x;x;x],[y;y;y;y],[z;z;z;z]}; %last line needed by all three

set_cont = {[1,2,3],[1,2,3],[1,2,3],1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3};

%clear old response functions and perm dipole shifts
clear R1 R2 R3 R4 R5 R6 pdm_shift pdm
if ~exist('params_for_wrapper','var')
    params_for_wrapper = [];
end
if ~exist('t_verbose','var')
    t_verbose = false; %don't output time as default
end
    
if ~exist(params_for_wrapper,'file')==2  
error('requires variable "params_for_wrapper" to be passed')
 end
    %if this variable is present and it is a string corresponding to a file
    %it uses this.  Clear the name if this isn't desired and you wish to
    %enter them into the wrapper manually
    
max_tier = 2; Kappa = 3; num_realisations = 1; %default values if nothing loaded
    
 load(params_for_wrapper,'Temperature_in_Kelvin','system_parameters_file',...
     'beam_param_set','om1_rng','t2_range_fs','om3_rng','pol','bb',...
    'max_tier','Kappa','num_realisations','t1_coh','om2_coh','om3_coh')

%pol = [x;y;x;y ];

if ~exist('pol','var');  loop_all_pol = true; %assume loop over all required
elseif islogical('pol');  loop_all_pol = true; %assume loop over all required    
else;  loop_all_pol = false; %take just one
end

L3 = length(om3_rng);  L2 = length(t2_range_fs); 
if iscell(om1_rng); L1 = sum(cellfun(@length, om1_rng)); 
else; L1 =length(om1_rng); end

    calc_rp = true;  calc_nr = true;  calc_coh = false; %not saved 
    [t_scale, B,speed_unit]= inv_cm_unit_sys(Temperature_in_Kelvin);
    t2_range = t2_range_fs*t_scale/1e3;   %in inverse cm 
    
    max_tier_final = max_tier;

%% Pass these parameters to the function which calculates quantities 
%independent of static disorder

[coup_com_save,coup_acom_save,const_factor,QQ_topass,nn,fock_space_rep,...
    H_site,site_shift,g_cnst,lam_tot,V_n] = ...
    TwoD_with_HEOM_prerun(system_parameters_file,Temperature_in_Kelvin,...
    num_realisations,max_tier,Kappa,[]) ;

%{
if conv_w_gaussian
%subtract mean difference so average excitation is the same, just include
%the extra shift at the end.
site_shift = site_shift - repmat(sum(site_shift,2),1,length(sd_mat))/N; 
end
    %}
%also find pdm_shift is an N X N cell array with the
%shift of the excited state with sites j and k excited from the PDM
load(system_parameters_file,'pdm_shift');
pdm_shift{1,2} =0; pdm_shift{2,1} =0;
%[LASTMSG, LASTID] = lastwarn; %turn off warnings that this variable hasn't loaded
%warning('off',LASTID)
if ~exist('pdm_shift','var') 
   load(system_parameters_file,'pdm');
   if exist('pdm','var') %perm dipole moments known but not shifts to states
        load(system_parameters_file,'R');
        for j =1:N
            for kk = 1:N
pdm_shift{kk,j} = dot(pdm(kk,:),pdm(j,:))-3*dot(pdm(kk,:),R(j,:))...
              *dot(pdm(j,:),R(kk,:))/norm(R(j,:)-R(kk,:))^2;
pdm_shift{kk,j} = pdm_shift/norm(R(j,:)-R(kk,:))^3  ;   
            end
        end       %loop over states
   else
       pdm_shift=[]; %no perm dipole moment
   end
end

HEOM_params = {coup_com_save,coup_acom_save,const_factor,...
             QQ_topass,nn,max_tier_final};
clear coup_com_save coup_acom_save const_factor QQ_topass nn max_tier_final

Ham_params = {fock_space_rep,H_site,site_shift,pdm_shift};

comp_lg = {calc_rp,calc_nr,calc_coh};


%% Now calculate the response functions
sdlp_rng = 1:num_realisations; % in general one can do this multiple 
%times until convergence is achieved

N = length(H_site); %number of chromophore sites
N2 = N*(N-1)/2; %number of double excited states
%calculate the operator for every interaction with the beam
lg = exist('bb','var');
if loop_all_pol 
       lp_rng = 1:length( required_pol);
    %preallocate the 3 different LI components with a 5th dimension for
    %magnetic induced chiral comps if bb is given as true
    
    sz5 =1; if lg; if islogical(bb) && bb; sz5 =5; else lg=false; end; end;
    if calc_nr; R1_s = complex(zeros(L3,L2,L1,3,sz5)); R4_s = R1_s;   R6_s = R1_s; end;
    if calc_rp; R2_s = complex(zeros(L3,L2,L1,3,sz5)); R3_s = R2_s;   R5_s = R2_s; end;
    pol_set = required_pol; bb_set = required_pol; %set as same
else
    lp_rng = 1; pol_set = {pol};   bb_set = {bb};
end

%should be the same after a random rotation
RR = eye(3);
RR= rot_3D(rand(3,1)); %this isn't a uniformly random rotation 
%but will do comment to remove but seems to make no difference

if lg
load(system_parameters_file,'mu','R','mdip')
mm = mdip - pi*1i*cross(R.',mu.').'; 
%note this is missing a factor of 1/(wavelength) as this depends on the
%light  lambda omega = 2*pi*c = 1 in our units so 1/lambda = omega
% Only 1st and 4th int are def in freq space
mm = mm * mean(om3_rng); %scale with mean output frequency
end
mu = (RR*mu.').'; mm = (RR*mm.').';

for lpp=lp_rng
lpp %output so as to follow this, delete this if you want to avoid seeing this
pol = pol_set(:,:,lpp);

V1 = zeros(N,1); [V2,V3,V4]=deal(V1,V1,V1);
VV1 = zeros(N2,N); [VV2,VV3,VV4]=deal(VV1,VV1,VV1);
%first those with nonchiral origins

inc_lg = mu*pol'; %should be N by 4
%if any column is all zeros then there is no contribution
inc_lg2 = all(inc_lg == 0);
for j = 1:N

    V1 = V1 + V_n(2:N+1,1,j).*inc_lg(j,1);
    V2 = V2 + V_n(2:N+1,1,j).*inc_lg(j,2);
    V3 = V3 + V_n(2:N+1,1,j).*inc_lg(j,3);
    V4 = V4 + V_n(2:N+1,1,j).*inc_lg(j,4);
  
  %  VV1 = VV1 + V_n(N+2:end,2:N+1,j).*dot(mu(j,:),pol(1,:)); %not required
    VV2 = VV2 + V_n(N+2:end,2:N+1,j).*inc_lg(j,2);
    VV3 = VV3 + V_n(N+2:end,2:N+1,j).*inc_lg(j,3);
    VV4 = VV4 + V_n(N+2:end,2:N+1,j).*inc_lg(j,4);
    
end
%now the operators are scaler rather than tensor
V_eg = {V1,V2,V3,V4}; V_fe = {[],VV2,VV3,VV4};
if ~any(inc_lg2) %if any of the operators are all zero this can be skipped
[R1,R2,R3,R4,R5,R6,lin_spec,tmp_gg_t,tmp_ee_t]=twoD_with_HEOM_NOA...
        (sdlp_rng,HEOM_params,Ham_params,comp_lg,V_eg,V_fe...
        ,om1_rng,t2_range,om3_rng,t_verbose);

if loop_all_pol
    for lp2 = set_cont{lpp}
        if calc_nr
            R1_s(:,:,:,lp2,1) = R1_s(:,:,:,lp2,1)+R1; 
            R4_s(:,:,:,lp2,1) = R4_s(:,:,:,lp2,1)+R4; 
            R6_s(:,:,:,lp2,1) = R6_s(:,:,:,lp2,1)+R6; 
        end
        if calc_rp
            R2_s(:,:,:,lp2,1) = R2_s(:,:,:,lp2,1)+R2; 
            R3_s(:,:,:,lp2,1) = R3_s(:,:,:,lp2,1)+R3; 
            R5_s(:,:,:,lp2,1) = R5_s(:,:,:,lp2,1)+R5; 
        end
%R_s in principle contains ALL information about the response function
% of the isotropic system
    end
end
end
    if lg %also calculate contributions from the magnetic part

        bb2 = bb_set(:,:,lpp);

        for mag_lp = 1:4
            %check other nonchiral interactions are non zero
        if ~any(inc_lg2([1:mag_lp-1,mag_lp+1:4]))     
            
        Vm_eg = V_eg; Vm_fe = V_fe; Vm = zeros(N,1); VVm = zeros(N2,N); 
        inc_lg3 = mm*bb2(mag_lp,:)';
        if ~all(inc_lg3 == 0) %check at least one non zero
             for j=1:N
                Vm = Vm + V_n(2:N+1,1,j).*inc_lg3(j);
                VVm = VVm + V_n(N+2:end,2:N+1,j).*inc_lg3(j);
             end
             %replace the mag_lpth dipole interaction with magnetic int
        Vm_eg{mag_lp} = Vm;    Vm_fe{mag_lp} = VVm; 
        
        [R1m,R2m,R3m,R4m,R5m,R6m,lin_spec,tmp_gg_t,tmp_ee_t]=twoD_with_HEOM_NOA...
        (sdlp_rng,HEOM_params,Ham_params,comp_lg,Vm_eg,Vm_fe...
        ,om1_rng,t2_range,om3_rng,t_verbose);

if loop_all_pol
    for lp2 = set_cont{lpp}
        if calc_nr
            R1_s(:,:,:,lp2,1+mag_lp) = R1_s(:,:,:,lp2,1+mag_lp)+R1m; 
            R4_s(:,:,:,lp2,1+mag_lp) = R4_s(:,:,:,lp2,1+mag_lp)+R4m; 
            R6_s(:,:,:,lp2,1+mag_lp) = R6_s(:,:,:,lp2,1+mag_lp)+R6m; 
        end
        if calc_rp
            R2_s(:,:,:,lp2,1+mag_lp) = R2_s(:,:,:,lp2,1+mag_lp)+R2m; 
            R3_s(:,:,:,lp2,1+mag_lp) = R3_s(:,:,:,lp2,1+mag_lp)+R3m; 
            R5_s(:,:,:,lp2,1+mag_lp) = R5_s(:,:,:,lp2,1+mag_lp)+R5m; 
        end
%elements 2:5 in R_s 
    end
end
        end
        end
        end
    end

end

%To calculate from R_s = [A;B;C] dot product it with 1/30*
% [(p1.p2)(p3.p4),(p1.p3)(p2.p4),(p2.p3)(p1.p4)]*[4,-1,-1;-1,4,-1;4,-1,-1]
    
%% Convolve with Gaussian if required
%doesn't yet work
%{
% Mean frequency shift Delta om_eg = sum_k Delta om_kg is also distributed
% Gaussianly, and is simple to include
N = length(H_site);
sd_var = sqrt(dot(sd_mat,sd_mat))/N; %sigma of normal dist of mean freq shift

%check there is actually static D otherwise no point
if sd_var == 0;  conv_w_gaussian = false;    end

if conv_w_gaussian 
load(system_parameters_file,'sd_mat');

zz = linspace(-10*sd_var,10*sd_var,4000);

Conv_dist = exp(-zz.^2/2/sd_var^2)/sqrt(2*pi*sd_var^2);
%Need to calculate 
%S'(om_3,t,om_1) = int_{-inf}^{inf} dw S(om_3+w,t,om_1+w)*P(w)
%%
%for frequency space rephasing signal
Conv_dist = full(sparse(1:length(Conv_dist),1:length(Conv_dist),Conv_dist));%diag(Conv_dist);

R1 = convn(R1,Conv_dist,'same');
R2 = convn(R2,Conv_dist,'same');
R3 = convn(R3,Conv_dist,'same');
R4 = convn(R4,Conv_dist,'same');
R5 = convn(R5,Conv_dist,'same');
R6 = convn(R6,Conv_dist,'same');

%% Inverse FFT for each value of t2, multiply by xi(t1 +/- t3) fft back

%can pass 
if iscell(om1_rng); t1_rng = om1_rng{2};om1_rng = om1_rng{1};  end;
if iscell(om3_rng); t3_rng = om3_rng{2};om3_rng = om3_rng{1};  end;

test1 = diff(om1_rng,2); test2 = diff(om3_rng,2); 
if any(abs(test1)>eps(max(om1_rng))) ||  any(abs(test2)>eps(max(om3_rng)))
warning('om1_rng and om3_rng are not equally spaced, this will cause problems')
end

N1 = length(om1_rng); N3 = length(om3_rng);
dt1 = 2*pi/(om1_rng(end)-om1_rng(1)); t1 = 0:dt1:(N1-1)*dt1;
dt3 = 2*pi/(om3_rng(end)-om3_rng(1)); t3 = 0:dt3:(N3-1)*dt3;

for lp =1:6 %signals
    for lp2 = 1:(size(av_4,5)+size(av_5,6)) %beam configurations
    Rbroad = zeros(length(om1_rng),length(t2_rng),length(om3_rng));
    eval(strcat('Rtmp=R',num2str(lp),'(:,:,:,',num2str(lp2),');')) 

        for t2lp = 1:length(t2_range)
            Rtmp2 = squeeze(Rtmp(1:length(om3_rng),t2lp,1:length(om1_rng)));
            Rtmp2 = ifft2(Rtmp2); %take inverse fft
            %correct for the fact that om1/3 don't start from zero
            expfact = exp(1i*(t1+t3)) 
            'unfinished'
        end
        
eval(strcat('R',num2str(lp),'(:,:,:,',num2str(lp2),')=Rbroad;'))   %replace with the correct thing        
    end
end

%% Shitty method
%{
om3m = min(om3_rng); om3M = max(om3_rng);
om1m = min(om1_rng); om1M = max(om1_rng);

PP = @(w) exp(-w.^2/2/sd_var^2)/sqrt(2*pi*sd_var^2);
rng1 = 1:length(t2_range_fs);

for lp =1:6
eval(strcat('Rbroad=R',num2str(lp),';')) 
if~isempty(Rbroad) %don't bother if I haven't set this
    Rtmp = Rbroad;  rng2 = 1:size(Rtmp,4);
    for lp1 = 1:length(om1_rng); om1 = om1_rng(lp);
        for lp3 = 1:length(om3_rng);  om3 = om1_rng(lp);

wmin = max(om3-om3m,om1-om1m); wmax =  min(om3M-om3,om1M-om1);
%range of w which I can use given the range of om3 and om1

norm_fct = integral(PP,wmin,wmax);
%effective normalisation due to trunction, should be ~ 1 in most cases and
%greater than a half
if length(rng2)>1
FF = @(w) PP(w).*interpn(om3_rng,rng1,om1_rng,rng2,Rtmp,om3+w,rng1,om1+w,rng2);
else
FF = @(x) PP(x).*interpn(om3_rng,rng1,om1_rng,Rtmp,om3+x,rng1,om1+x);    
end
Rbroad(lp3,:,lp1,:) = integral(FF,wmin,wmax,'RelTol',1e-6,'ArrayValued',true)/norm_fct; %#ok<SAGROW>
        end
    end
eval(strcat('R',num2str(lp),'=Rbroad;'))   %replace with the correct thing
end
end
clear Rtmp Rbroad
%}
end
%}