% This code calculates 2D spectroscopy using the HEOM in the exciton basis
% for an excitonic system with a spectral density of underdamped and
% overdamped brownian oscillator modes on each side.

%clear old response functions and perm dipole shifts
clear R1 R2 R3 R4 R5 R6 pdm_shift pdm
if ~exist('params_for_wrapper','var')
    params_for_wrapper = [];
end

if exist(params_for_wrapper,'file')==2  
    %if this variable is present and it is a string corresponding to a file
    %it uses this.  Clear the name if this isn't desired and you wish to
    %enter them into the wrapper manually
num_realisations = 1; %default values if nothing loaded
  pop_only = true; %only consider exciton-exciton transitions
 load(params_for_wrapper,'Temperature_in_Kelvin','pop_only','system_parameters_file',...
     'beam_param_set','om1_rng','t2_range_fs','om3_rng','t1_coh','om2_coh','om3_coh',...
    'num_realisations','numvib')

    calc_rp = true;  calc_nr = true;  calc_coh = false; %not saved 
    [t_scale, B,speed_unit]= inv_cm_unit_sys(Temperature_in_Kelvin);
    t2_range = t2_range_fs*t_scale/1e3;   %in inverse cm 

else %enter them into this format into the wrapper
chooser = 2; %choose which system stuff to use, this is just a lazy setup 

%% System related parameters
%enter the name of the file where the parameters are stored
if chooser==0
system_parameters_file = 'Parameters/dimer_system_params.mat';
numvib{1} = 3; numvib{2} = 3; 
Temperature_in_Kelvin = 300; %self explainatory
elseif chooser  == 1
system_parameters_file = 'Parameters/Adv_HEOM_for_Efficient_eval_params.mat'; 
numvib{1} = 3; numvib{2} = 3; 
Temperature_in_Kelvin = 77; %self explainatory    
else %to reproduce for alexandra etc etc
system_parameters_file = 'Parameters/plenio_paper_parameters.mat'; 
numvib{1} = 3; numvib{2} = 3; 
Temperature_in_Kelvin = 77; %self explainatory
end
%calculate scaling from ps->inv cm, thermo beta in cm^(-1) 
[t_scale, B,speed_unit]= inv_cm_unit_sys(Temperature_in_Kelvin);
 %multiply time in ps by t_scale to get it in inverse cm
 %t_in_ps * t_scale == t_in_inverse_cm
 
num_realisations = 1; %realisations of static disorder

%% Beam params, wavevector and polarization

k1 = [1/sqrt(2),0,1/sqrt(2)]; pol_1 = [0,1,0];  mag_1 = cross(k1,pol_1);
k2 = [1/sqrt(2),0,-1/sqrt(2)]; pol_2 = [0,1,0]; mag_2 = cross(k2,pol_2); 
k3 = [0,0,1]; pol_3 = [0,1,0]; mag_3 = cross(k3,pol_2); 
%also specify the signal used for Heterodyne detection (will be conjugated)
 pol_HET = [0,1,0];  pol_HET2 = [1,0,0]; 

beam_param_set = cat(3,[pol_1;pol_2;pol_3;pol_HET;k1;k2;k3],...
                 [pol_1;pol_2;pol_3;pol_HET2;k1;k2;k3]);
 %put all parameters into a single set 
             
calc_rp = true;  calc_nr = true;  calc_coh = false; 

beam_param_set = [1,0,0;1,0,0;1,0,0;1,0,0;0,0,1;0,0,1;0,0,1];

%% Choose the range of frequencies to tabulate for RP/NRP
%Calculate R(om_3,t_2,om_1) = iint_dt_3 dt_1 R(t_3,t_2,t_1)
if chooser ==0
 om1_rng = linspace(10^7./(850),10^7./(730),100); %range in cm^(-1)
t2_range_fs = (0:2:60);  %delay time t_2 in fs
om3_rng = linspace(10^7./(920),10^7./(680),1000); 
t2_range = t2_range_fs*t_scale/1e3;   %in inverse cm
elseif chooser == 1

om1_rng = linspace(10000-800,10000+800,100); %range in cm^(-1)
t2_range_fs = [0,100,200,300,400,500];  %delay time t_2 in fs
om3_rng = linspace(10000-800,10000+800,1000); 
t2_range = t2_range_fs*t_scale/1e3;    %in inverse cm
else
    
om1_rng = linspace(1.21e4,1.27e4,100); %range in cm^(-1)

t2_range = linspace(0,0.1586,100);
t2_range_fs = t2_range/t_scale*1e3;  %delay time t_2 in fs
om3_rng = linspace(1.21e4,1.27e4,400) ;
    
end
%{
% dom = om1_rng(2)-om1_rng(1);
% om_end = om1_rng(end)-om1_rng(1);
% t_ef = 2*pi*(0:1/om_end:1/dom);
%take long well spaced time intervals when calculating the coherences,
%possible as I am using FFT
t1_range = (0:1e-4:2) /t_scale; 
t3_range = (0:1e-4:2) /t_scale;
%}
%% If required choose range for the coherence pathway

t1_coh = 40; %pass single number to get the program to choose a range 
%based on markovian decay with this number of points
om2_coh= linspace(10^7./(1000),10^7./(700),100); %range 
om3_coh= linspace(10^7./(1000),10^7./(700),100); %range 

end
%% Pass these parameters to the function which calculates quantities 
%independent of static disorder

[fock_space_rep,H_el,H_vib,H_e,H_f,site_shift,...
    V_n,V2_n,av_2,av_3,av2_3,av_4,av_5,av2_5,Drude_modes,BO_modes] = ...
    TwoD_with_Redfield_prerun...
    (system_parameters_file,num_realisations,beam_param_set,numvib) ;

%% Now calculate the response functions

sdlp_rng = 1:num_realisations; % in general one can do this multiple 
%times until convergence is achieved
if size(Drude_modes{1},2)==2
opts = 'Redfield' ;   
else
opts = 'Pure_dephasing';
end

[R1,R2,R3,R4,R5,R6,lin_spec,tmp_gg_t,tmp_ee_t,R_red_op]=twoD_with_Redfield_main2(...
        sdlp_rng,fock_space_rep,H_el,H_vib,H_e,H_f,site_shift,B...
       ,Drude_modes,BO_modes,V_n,av_2,av_3,av2_3,av_4,av_5,av2_5,  ...
        calc_rp,calc_nr,om1_rng,t2_range,om3_rng,opts,true);
   
%{
load(system_parameters_file,'mu');
[R1,R2,R3,R4,R5,R6,lin_spec,tmp_gg_t,tmp_ee_t]=twoD_with_Redfield_main(...
        sdlp_rng,fock_space_rep,H_el,H_vib,H_e,H_f,site_shift,B...
       ,Drude_modes,BO_modes,V_n,V2_n,mu,{[],[]},beam_param_set,[],[],[],  ...
        calc_rp,calc_nr,om1_rng,t2_range,om3_rng,opts,true);     
        %}
     