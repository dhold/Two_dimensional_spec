function [Sig_save,pt] = pp_extract(save_file,cmp, pnts1,pnts2...
                           ,pulse_sd_cm,t_scale,tau_rng)

%gets pump probe data from response function S(om_s,tau,om_1) saved in a 
%file save_file, this assumes that the probe is fast and frequency resolved
%detection is applied.  If tau_range
%pnts1 %carrier frequency of omega_1 to select
%pnts2 %frequency of omega_sig to select in cm^{-1}
%pulse_sd_cm = sigma in E(t) propto exp(-(t-tau)^2/2/sigma^2) 
%t_scale scales from fs to inverse cm, 
%tau_rng is the range of separations you wish to consider, passing this
%argument will NOT ignore time dependence over t2 during the interaction

fle = matfile(save_file);
om3_rng = fle.om3_rng;  om1_rng = fle.om1_rng; 
t2_range = fle.t2_range_fs/t_scale;  %t_scale converts to inverse cm

if ~exist('tau_rng','var')
    % Take integral over range corresponding to many frequencies
%ignoring the time dependence over t2
GSB = imag(fle.R3(:,:,:,cmp)+fle.R4(:,:,:,cmp));
SE = imag(fle.R2(:,:,:,cmp)+fle.R1(:,:,:,cmp));
ESA = imag(fle.R5(:,:,:,cmp)+fle.R6(:,:,:,cmp));

for lp =1:length(pnts1)
    %pulse frequency envelope with the carrier frequency
    E_w_sq = exp(-pulse_sd_cm^2*(pnts1(lp)-om1_rng).^2)*sqrt(pulse_sd_cm^2/pi);
    E_w_sq = repmat(reshape(E_w_sq,length(E_w_sq),1),[1,length(pnts2),length(t2_range_fs)]);
    tmp = interp1(om3_rng,GSB+SE+ESA,pnts2);
Sig_save{lp} = squeeze(trapz(om1_rng,E_w_sq.*permute(tmp,[3,1,2])));
end
else
% Take integral NOT ignoring the time dependence over t2
%this is mucchhhh slower and requires loading rephasing and nonrephasing
%separately
%find nearest points to these, I cba with interpolation honestly
pt = zeros(length(pnts2),1);
for lp3 =1:length(pnts2)
[~,pt(lp3)] = min(abs(pnts2(lp3)-om3_rng)); 
end
Srpfull = fle.R3(pt,:,:,cmp)+ fle.R2(pt,:,:,cmp)+ fle.R5(pt,:,:,cmp);
Snrfull = fle.R4(pt,:,:,cmp)+ fle.R1(pt,:,:,cmp)+ fle.R6(pt,:,:,cmp);

[ww,tt] = meshgrid(om1_rng,t2_range); %meshgrid points
Sig_save = zeros(length(tau_rng),length(pnts1),length(pnts2));
%find nearest points to these, I cba with interpolation honestly
for lp3 =1:length(pnts2)

Srp = squeeze(Srpfull(lp3,:,:)); Snr = squeeze(Snrfull(lp3,:,:));
%take real and imag because of the fft

for lp =1:length(pnts1)
    ww2 = (pnts1(lp)-ww);
    E_w = exp(-pulse_sd_cm^2*ww2.^2/2)*(pulse_sd_cm^2/pi)^(1/4);
    for tau_lp = 1:length(tau_rng)
        t_sep = tau_rng(tau_lp);  tt2 = t_sep-tt;
    %pulse frequency envelope with the carrier frequency
    
    E_t = exp(-tt2.^2/2/pulse_sd_cm^2)/(pulse_sd_cm^2*pi)^(1/4);
  %  E_t = E_t.*heaviside(t_sep-tt);

    tmp = trapz(t2_range,E_w.*E_t.*(exp(1i.*ww2.*tt2).*Snr + exp(-1i.*ww2.*tt2).*Srp));
    tmp = trapz(om1_rng,tmp);
   % tmp2 = simp2D((exp(1i.*ww.*tt).*Snr + exp(-1i.*ww.*tt).*Srp)...
    %               .* E_w.*E_t ,om0,ome,t0,te);  
     Sig_save(tau_lp,lp,lp3) = tmp;
    end
end
end
end

%{
fle = matfile('saved_data\PPcd_with_HEOM_no_mode.mat');
om1_rng = fle.om1_rng;  om1_mid = mean(om1_rng); dom_1 = om1_rng(2)-om1_rng(1);
om3_rng = fle.om3_rng;  om3_mid = mean(om3_rng); dom_3 = om3_rng(2)-om3_rng(1);
[t_scale]= inv_cm_unit_sys(0);

%R1 = squeeze(fle.R1(:,:,:,cmp));   R2 = squeeze(fle.R2(:,:,:,cmp));

%% Load Rephasing and non rephasing
tau = 150; tau = tau*t_scale/1000;
t2_range = t2_range_fs*t_scale/1000;
om_pump_rng = 1e4*[1.2255,1.2400,1.2450,1.2500]; %pump values to consider

SE = zeros(length(om3_rng),length(t2_range),length(om_pump_rng));

N1 = 101; N2 = 101; Nsd = 7;

tprime_rng = tau*linspace(-Nsd,Nsd,N1); %time 
om_rng = linspace(-Nsd,Nsd,N2)/tau; %freq, take NSD each way

% Set up the integration step sizes in the x and y directions
hx = 2*Nsd*tau/N1;
hy = 2*Nsd/tau/N2;

% define grid vectors

[tt,ww] = meshgrid( tprime_rng ,om_rng );

EE = @(t) exp(-t.^2/tau^2/2)/sqrt(sqrt(pi)*tau);
EE2 = @(w) exp(-w.^2*tau^2/2)/sqrt(sqrt(pi)/tau);

Efact1 = EE(tt).*EE2(ww).*exp(1i*tt.*ww);
Efact2 = EE(tt).*EE2(ww).*exp(-1i*tt.*ww);

for tau_lp = 1:10%length(t2_range)
    t_sep = t2_range(tau_lp); 
   % lg = t_sep - tprime_rng >= 0  & t_sep - tprime_rng <=max(t2_range) ;
    for om3_lp = 1:length(om3_rng)
        for om1_lp = 1:length(om_pump_rng)
            om_u = om_pump_rng(om1_lp)-om1_mid;
            
            %load and interpolate to the correct range
R1tmp = squeeze(R1(om3_lp,:,:));  
R1tmp = interpn(t2_range,om1_rng-om1_mid,R1tmp,t_sep+tt,om_u-ww,'linear',0);
R2tmp = squeeze(R2(om3_lp,:,:));
R2tmp = interpn(t2_range,om1_rng-om1_mid,R2tmp,t_sep+tt,om_u+ww,'linear',0);

U = Efact2.*R2tmp + Efact1.*R2tmp;

SE(om3_lp,tau_lp,om1_lp) = imag(simp2D(U,N1-1,N2-1,hx,hy)); %integrate over the range using simpsons rule

        end
    end
end
%}


