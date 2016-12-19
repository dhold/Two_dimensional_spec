%%
% om1_rng = linspace(1.21e4,1.27e4,100); %range in cm^(-1)
% om3_rng = linspace(1.21e4,1.27e4,300) ;
t2_range_fs = linspace(0,2000,210); 
CMRmapsmall=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
CMRmap=interp1(1:9,CMRmapsmall,1:0.5:9);


%fle = matfile('2Dcd_with_HEOM_no_mode.mat')
%fle = matfile('2Dcd_with_HEOM_with_mode.mat');
%fle = matfile('2Dcd_with_HEOM_weaker_mode2.mat');
%fle = matfile(file_name);
%fle = matfile('saved_data/PPcd_with_HEOM_with_mode');
fle = matfile('saved_data/2Dcd_HEOM_with_mode.mat');
%fle = matfile('saved_data/2Dcd_with_RF_and_mode.mat');
%fle = matfile('saved_data/PPcd_with_HEOM_pop_only_no_mode.mat');
%fle = matfile('saved_data/PPcd_with_HEOM_fastpump_no_mode_.mat');
om1_rng = fle.om1_rng; if iscell(om1_rng); om1_rng = om1_rng{1}; end
rnj = 1:length(om1_rng);
om3_rng = fle.om3_rng;

scale_fct2 = max(max(max(abs(fle.R3(:,1,rnj,1)+fle.R4(:,1,rnj,1)))));

%%
bpset = fle.beam_param_set; ps = 1; t2L = length(t2_range_fs);
cmp = size(bpset{1},3)+ps; tmp = bpset{2}(:,:,ps);
cmp2 = size(bpset{1},3)+2; %component to subtract from primary, set to 0 for none
fct = +1;
%{
figure; contourf(om1_rng/1e4,om3_rng/1e4,imag(squeeze(...
fle.R2(:,t2L,:,5)+fle.R2(:,t2L,:,4)+fle.R5(:,t2L,:,5)+fle.R5(:,t2L,:,4)))/scale_fct2); 
figure; contourf(om1_rng/1e4,om3_rng/1e4,imag(squeeze(...
fle.R2(:,t2L,:,4)+fle.R2(:,t2L,:,7)+fle.R5(:,t2L,:,4)+fle.R5(:,t2L,:,7)))/scale_fct2);

figure; contourf(om1_rng/1e4,om3_rng/1e4,imag(squeeze...
    (fle.R2(:,t2L,:,5)+fle.R2(:,t2L,:,8)+fle.R5(:,t2L,:,5)+fle.R5(:,t2L,:,8)))/scale_fct2); 
%%
figure; contourf(om1_rng/1e4,om3_rng/1e4,imag(squeeze(...
fle.R2(:,t2L,:,6)-fle.R2(:,t2L,:,9)))/scale_fct2); 
figure; contourf(om1_rng/1e4,om3_rng/1e4,imag(squeeze(...
fle.R2(:,t2L,:,7)-fle.R2(:,t2L,:,8)))/scale_fct2); 
figure; contourf(om1_rng/1e4,om3_rng/1e4,imag(squeeze(...
fle.R2(:,t2L,:,3+8)-fle.R2(:,t2L,:,3+9)+0*fle.R5(:,t2L,:,3+8)-0*fle.R5(:,t2L,:,3+9)))/scale_fct2); 
%}
%% Pumpprobe
GSB = imag(fle.R3(:,:,rnj,cmp)+fle.R4(:,:,rnj,cmp))/scale_fct2; 
SE = imag(fle.R2(:,:,rnj,cmp)+fle.R1(:,:,rnj,cmp))/scale_fct2; 
ESA = imag(fle.R5(:,:,rnj,cmp)+fle.R6(:,:,rnj,cmp))/scale_fct2;
sig_type= 3;
data_set = 'pump probe'
%% Rephasing
GSB = imag(fle.R3(:,:,rnj,cmp))/scale_fct2; 
SE = imag(fle.R2(:,:,rnj,cmp))/scale_fct2; 
ESA = imag(fle.R5(:,:,rnj,cmp))/scale_fct2;

if cmp2 ~= 0
GSB =GSB+fct*imag(fle.R3(:,:,rnj,cmp2))/scale_fct2; 
SE = SE+fct*imag(fle.R2(:,:,rnj,cmp2))/scale_fct2; 
ESA = ESA+fct*imag(fle.R5(:,:,rnj,cmp2))/scale_fct2;
end
sig_type= 1;
data_set = 'rephasing'
fle.chk(ps,1) %check this is actually a valid configuration
%% NonRephasing
GSB = imag(fle.R4(:,:,rnj,cmp))/scale_fct2; 
SE = imag(fle.R1(:,:,rnj,cmp))/scale_fct2; 
ESA = imag(fle.R6(:,:,rnj,cmp))/scale_fct2;
sig_type= 2;
data_set = 'nonrephasing'
fle.chk(ps,2)
%% Calc max and min for colour bars

min1= min(GSB(:)); max1=max(GSB(:));
min2= min(SE(:)); max2=max(SE(:));
min3= min(ESA(:)); max3=max(ESA(:));

%%
figure
title(strcat('t=',num2str(t2_range_fs),'fs'))
contourf(om1_rng/1e4,om3_rng/1e4,squeeze(0*GSB(:,1,:)+GSB(:,end,:)),size(CMRmap,1))
 colormap(CMRmap)
%shading interp
colorbar
%caxis([min1,max1])
xlabel('\omega_1  (10^{4} cm^{-1})')
ylabel('\omega_s  (10^{4} cm^{-1})')
%% Create multi figure
min4 = min(min2,min3); max4 = max(max2,max3);
for pnt = 1:10:71

t2_range_fs(pnt)
om1lg = om1_rng>12100 & om1_rng<12650;
om3lg = om3_rng>12100 & om3_rng<12650;

figure;   set(gcf, 'renderer', 'zbuffer');    subplot(2,1,1); 
contourf(om1_rng(om1lg)/1e4,om3_rng(om3lg)/1e4,squeeze(SE(om3lg,pnt,om1lg)),size(CMRmap,1))
 colormap(CMRmap);   colorbar
%caxis([min4,max4])
ylabel('\omega_s  (10^{4} cm^{-1})','FontSize',14)

subplot(2,1,2);  
contourf(om1_rng(om1lg)/1e4,om3_rng(om3lg)/1e4,squeeze(ESA(om3lg,pnt,om1lg)),size(CMRmap,1))
colorbar
%caxis([min4,max4])
 colormap(CMRmap) %shading interp
xlabel('\omega_1  (10^{4} cm^{-1})','FontSize',14)
ylabel('\omega_s  (10^{4} cm^{-1})','FontSize',14)

%[a,b] = find(tmp'); 
fig_name = strcat('ES_NR_xxyx_zzz-yyz',num2str(t2_range_fs(pnt),3),'fs.png');
set(gcf, 'Color', 'w');
%export_fig(fig_name,'-p0.01','-native') %pad slightly to avoid over cropping  
end
%% Calculate scaling
tmp_max = t2_range_fs*0;
for lp =1:length(t2_range_fs)
tmp_max(lp) = max(max(abs(squeeze(SE(:,lp,:)+ESA(:,lp,:)))));
end
scale_plot = true;
f = ezfit(tmp_max,'exp'); %fit exp decay to this
%also account for offsets using the first parameters as an initial guess.  
f2 = ezfit(tmp_max,'a*exp(b*x)+c',[f.m(1),f.m(2),tmp_max(end)-f.m(1)*exp(f.m(2)*f.x(end))]);
if scale_plot
approx_sig_decay = f2.m(1)*exp(f2.m(2)*f.x)+f2.m(3);
min4= min(SE(:)+ESA(:))/(f2.m(1)+f2.m(3)); max4 =max(SE(:)+ESA(:))/(f2.m(1)+f2.m(3));
else
approx_sig_decay = ones(size(t2_range_fs));
min4= min(SE(:)+ESA(:));     max4 =max(SE(:)+ESA(:));
end
%(data - c)/a
%% Create animated plot

numframes = floor(length(t2_range_fs)/2);
M(1:numframes) = struct('cdata',[],'colormap',[]);
figure

%GSB_change = GSB - repmat(GSB(:,1,:),[1,length(t2_range_fs),1]);
%min5 = min(GSB_change(:)); max5 = max(GSB_change(:));
for lp = 2:numframes
    
    contourf(om1_rng/1e4,om3_rng/1e4,squeeze(GSB(:,lp,:)+SE(:,lp,:)+ESA(:,lp,:)...
        -0*(GSB(:,1,:)+SE(:,1,:)+ESA(:,1,:)))./approx_sig_decay(lp),size(CMRmap,1))
    caxis([min4,max4])
% contourf(om1_rng/1e4,om3_rng/1e4,squeeze(GSB_change(:,lp,:)))
%  caxis([min5,max5])
   colormap(CMRmap)

colorbar
title(strcat('\tau=',num2str(t2_range_fs(lp)),'fs'))


M(lp) = getframe;
end
%
%% alternative
writerObj = VideoWriter('yxx+xyx_zzz_RP.avi');

open(writerObj);
%Generate initial data and set axes and figure properties.

figure
contourf(om1_rng/1e4,om3_rng/1e4,...
    squeeze(GSB(:,1,:)+SE(:,1,:)+ESA(:,1,:)),size(CMRmap,1));
axis tight
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');
numframes = floor(length(t2_range_fs)/4);

for lp = 2:numframes
    
    contourf(om1_rng/1e4,om3_rng/1e4,squeeze(GSB(:,lp,:)+...
        SE(:,lp,:)+ESA(:,lp,:))./approx_sig_decay(lp),size(CMRmap,1));
    caxis([min4,max4])
    colormap(CMRmap)
    colorbar
title(strcat('\tau=',num2str(t2_range_fs(lp)),'fs'))
for lp2=1:4
   frame = getframe;
   writeVideo(writerObj,frame);
end
end
%writerObj.FrameRate = 10;
close(writerObj);

%********** Next things actually take slices ****
%% Plot slice at a value of omega_1 

pnts1 = 12278; %frequency of omega_1 to select
pnts2 = [12278,12475];  %frequency of omega_sig to select
Sig = squeeze(interp1(om1_rng,permute(GSB + SE + ESA,[3,1,2]),pnts1));
Sig2 = interp1(om3_rng,Sig,pnts2);
figure; plot(t2_range_fs,Sig2)
%contourf(om3_rng/1e4,t2_range_fs,Sig,size(CMRmap,1))

%% Take integral over range corresponding to many frequencies
%ignoring the time dependence over t2
pulse_sd_fs = 150;  pulse_sd_cm = pulse_sd_fs/1e3*t_scale;
pnts1 = [12278,12475]; %carrier frequency of omega_1 to select
pnts2 = [12278,12475];  %frequency of omega_sig to select
clear Sig_save
for lp =1:length(pnts1)
    %pulse frequency envelope with the carrier frequency
    E_w_sq = exp(-pulse_sd_cm^2*(pnts1(lp)-om1_rng).^2)*sqrt(pulse_sd_cm^2/pi);
    E_w_sq = repmat(reshape(E_w_sq,length(E_w_sq),1),[1,length(pnts2),length(t2_range_fs)]);
    tmp = interp1(om3_rng,GSB+SE+ESA,pnts2);
Sig_save{lp} = squeeze(trapz(om1_rng,E_w_sq.*permute(tmp,[3,1,2])));
end
figure; plot(t2_range,Sig_save{1})
%% Take integral NOT ignoring the time dependence over t2
%this is mucchhhh slower and requires loading rephasing and nonrephasing
%separately

%fle2 = matfile('2Dcd_with_HEOM_no_mode.mat'); 
fle2 = matfile('saved_data/PPcd_with_HEOM_with_mode');
%fle2 = matfile('2Dcd_with_HEOM_weaker_mode2.mat');
tmp = fle2.beam_param_set;
cmp = size(tmp{1},3)+1; %select component corresponding to a valid PP setup e.g. xxyx zzz
om1_rng = fle2.om1_rng; om0 = om1_rng(1); ome = om1_rng(end);
t0 = t2_range(1); te = t2_range(end);
om3_rng = fle2.om3_rng;

pulse_sd_fs = 100;  pulse_sd_cm = pulse_sd_fs/1e3*t_scale;
pnts1 = [12278,12475]; %carrier frequency of omega_1 to select
tau_rng = linspace(0,1200,200)/1e3*t_scale; %range of separations
pnts2 = [12278,12475];  %frequency of omega_sig to select
%pnts1 = [12278,12460,12495]; pnts2 = [12270,12495];


[ww,tt] = meshgrid(om1_rng,t2_range); %meshgrid points
Sig_save = zeros(length(tau_rng),length(pnts1),length(pnts2));
%find nearest points to these, I cba with interpolation honestly
for lp3 =1:length(pnts2)
[~,pt] = min(abs(pnts2(lp3)-om3_rng)); 
%Srp = squeeze(fle2.R3(pt,:,:,cmp)+ fle2.R2(pt,:,:,cmp)+ fle2.R5(pt,:,:,cmp));
%Snr = squeeze(fle2.R4(pt,:,:,cmp)+ fle2.R1(pt,:,:,cmp)+ fle2.R6(pt,:,:,cmp));
%Srp = Srp-squeeze(fle2.R3(pt,:,:,cmp+1)+ fle2.R2(pt,:,:,cmp+1)+ fle2.R5(pt,:,:,cmp+1));
%Snr = Snr-squeeze(fle2.R4(pt,:,:,cmp+1)+ fle2.R1(pt,:,:,cmp+1)+ fle2.R6(pt,:,:,cmp+1));
Srp = squeeze(fle2.R3(pt,:,:,cmp)-fle2.R3(pt,:,:,cmp+1));
Snr = squeeze(fle2.R4(pt,:,:,cmp)-fle2.R4(pt,:,:,cmp+1));
%take imag and imag because of the fft

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
    % Sig_save2(tau_lp,lp,lp3) = tmp2;
    end
end
end

figure;plot0 = plot(tau_rng/t_scale*1e3,imag(Sig_save(:,:))/scale_fct2);
xlabel('t_2 (fs)','FontSize',16); ylabel('PP signal, a.u.','FontSize',14);
    for lp2 = 1:length(pnts2)
        for lp1 = 1:length(pnts1)
        dispnm = strcat('\omega_1 =',num2str(pnts1(lp1)),', \omega_2 =',num2str(pnts2(lp2)));
        set(plot0(lp1+(lp2-1)*(length(pnts1))),'DisplayName',dispnm)
        end
    end
legend('show')   
%%
    fig_name = strcat('xxyx_zzz-yyz_RP_with_mode_0fs.png');
set(gcf, 'Color', 'w');
export_fig(fig_name,'-p0.01','-painters')
%% Compute prony analysis at a point

 N_order = 9; lg = t2_range_fs < 1000; %choose range
t_red= t2_range_fs(lg ); dt2 = (t_red(2)- t_red(1)); fs = 1/dt2;
a_list = zeros(2,N_order); spoles= a_list; 
tau_list = spoles; omega_list = spoles;
 Temperature_in_Kelvin = 77;   
[t_scale, B, ~]= inv_cm_unit_sys(Temperature_in_Kelvin);
for lpval = 1:2
    if lpval ==1
       pnt = [1.2255,1.2492]*10^4; %pnt = [1.2492,1.2255]*10^4;% pnt = [1.2255,1.2492]*10^4; %pnt = [1.2492,1.2492]*10^4; %
    imp_resp = interpn(om3_rng,t_red,om1_rng,SE(:,lg,:),pnt(2),t_red,pnt(1));
    else
       pnt = [1.2255,1.2255]*10^4;%pnt = [1.2255,1.2492]*10^4;%pnt = [1.2492,1.2255]*10^4; %pnt = [1.2255,1.2255]*10^4;% 
    imp_resp = interpn(om3_rng,t_red,om1_rng,ESA(:,lg,:),pnt(2),t_red,pnt(1));
    end
[num,den] = prony(imp_resp,N_order,N_order);
[r,p,k]= residuez(num,den); %find residues

a_list(lpval,:)= r(:); %amplitude list
spoles(lpval,:)= log(p(:))*fs; %poles
tau_list(lpval,:)= 1./imag(spoles(lpval,:));%decay factor
omega_list(lpval,:) = imag(spoles(lpval,:));
%signals are imag so will be made of cosines

figure
 curve = zeros(N_order,length(t_red));
    for lp = 1:N_order
        if omega_list(lpval,lp) > eps
            fct = 2;
        elseif omega_list(lpval,lp) < -eps
            fct = NaN;
        else
            fct = 1;
        end
        curve(lp,:) = fct*a_list(lpval,lp).*exp(spoles(lpval,lp)*t_red);    
    end
    plot(t_red,imag(curve))
    xlabel('t_2 (fs)','FontSize',16); ylabel('ES signal, a.u.','FontSize',14);
    %for lp =1:N_order
    %    set('DisplayName',num2str(omega_list(lp)*1e3/t_scale))
    %end
end

%% compare signals at t=whatever
%fle = matfile(file_name);
fle = matfile('saved_data/2Dcd_with_HEOM_no_mode.mat'); %seems to be a
%sign error with this verion...
%fle = matfile('saved_data/2Dcd_HEOM_with_mode.mat');
tmp = fle.beam_param_set; pnt =1;

om1_rng = fle.om1_rng; if iscell(om1_rng); om1_rng = om1_rng{1}; end
rnj = 1:length(om1_rng);   om3_rng = fle.om3_rng;

om1lg = om1_rng>12100 & om1_rng<12750;
om3lg = om3_rng>12000 & om3_rng<12750;
tmp1 = 1:length(om1lg); tmp1 = tmp1(om1lg);
tmp3 = 1:length(om3lg); tmp3 = tmp3(om3lg);
sig_type = 1; scale_fct2 =1e-8;

tot_sig_set = [3,6,-1;1,2,1];

for cmp_lp = 1:2%size(tot_sig_set,1)
    
    cmp1 = size(tmp{1},3) +tot_sig_set(cmp_lp,1);
    cmp2 = size(tmp{1},3) +tot_sig_set(cmp_lp,2); 
    cnfg_s{cmp_lp} ={tmp{2}(:,:,tot_sig_set(cmp_lp,1)),tmp{2}(:,:,tot_sig_set(cmp_lp,2))};
    fct = tot_sig_set(cmp_lp,3);
    
if sig_type==1   
    
GSB2 = imag(fle.R3(tmp3,pnt,tmp1,cmp1)+fct*fle.R3(tmp3,pnt,tmp1,cmp2))/scale_fct2; 
SE2 = imag(fle.R2(tmp3,pnt,tmp1,cmp1)+fct*fle.R2(tmp3,pnt,tmp1,cmp2))/scale_fct2; 
ESA2 = imag(fle.R5(tmp3,pnt,tmp1,cmp1)+fct*fle.R5(tmp3,pnt,tmp1,cmp2))/scale_fct2;
%{
% GSB = imag(fle.R3(:,pnt,:,cmp))/scale_fct2; 
% SE = imag(fle.R2(:,pnt,:,cmp))/scale_fct2; 
% ESA = imag(fle.R5(:,pnt,:,cmp))/scale_fct2;
% GSB = GSB-imag(fle.R3(:,pnt,:,cmp+3))/scale_fct2; 
% SE = SE-imag(fle.R2(:,pnt,:,cmp+3))/scale_fct2; 
% ESA = ESA-imag(fle.R5(:,pnt,:,cmp+3))/scale_fct2;

%optical rotational part
% GSB = imag(fle.R3(:,pnt,:,cmp))/scale_fct2; 
% SE = imag(fle.R2(:,pnt,:,cmp))/scale_fct2; 
% ESA = imag(fle.R5(:,pnt,:,cmp))/scale_fct2; 

% GSB = imag(fle.R3(:,pnt,:,cmp))/scale_fct2; 
% SE = imag(fle.R2(:,pnt,:,cmp))/scale_fct2; 
% ESA = imag(fle.R5(:,pnt,:,cmp))/scale_fct2;
% GSB = GSB-imag(fle.R3(:,pnt,:,cmp+3))/scale_fct2; 
% SE = SE-imag(fle.R2(:,pnt,:,cmp+3))/scale_fct2; 
% ESA = ESA-imag(fle.R5(:,pnt,:,cmp+3))/scale_fct2;
%}
%'rephasing'
elseif sig_type==2 % NonRephasing
GSB2 = imag(fle.R4(tmp3,pnt,tmp1,cmp))/scale_fct2; 
SE2 = imag(fle.R1(tmp3,pnt,tmp1,cmp))/scale_fct2; 
ESA2 = imag(fle.R6(tmp3,pnt,tmp1,cmp))/scale_fct2;
else
GSB2 = imag(fle.R4(tmp3,pnt,tmp1,cmp)+fle.R3(tmp3,pnt,tmp1,cmp))/scale_fct2; 
SE2 = imag(fle.R1(tmp3,pnt,tmp1,cmp)+fle.R2(tmp3,pnt,tmp1,cmp))/scale_fct2; 
ESA2 = imag(fle.R6(tmp3,pnt,tmp1,cmp)+fle.R5(tmp3,pnt,tmp1,cmp))/scale_fct2;    
end
% Calc max and min for colour bars

min1= min(GSB2(:)); max1=max(GSB2(:));
min2= min(SE2(:)); max2=max(SE2(:));
min3= min(ESA2(:)); max3=max(ESA2(:));
% Create multi figure

figure
subplot(3,1,1); 
title(strcat('t=',num2str(t2_range_fs),'fs'))
contourf(om1_rng(om1lg)/1e4,om3_rng(om3lg)/1e4,squeeze(GSB2),size(CMRmap,1))
 colormap(CMRmap);  colorbar;  caxis([min1,max1]);

subplot(3,1,2); 
contourf(om1_rng(om1lg)/1e4,om3_rng(om3lg)/1e4,squeeze(SE2),size(CMRmap,1))
 colormap(CMRmap); colorbar;  caxis([min2,max2]);

subplot(3,1,3);  
contourf(om1_rng(om1lg)/1e4,om3_rng(om3lg)/1e4,squeeze(ESA2),size(CMRmap,1))
colormap(CMRmap); colorbar;  caxis([min3,max3]);  

xlabel('\omega_1  (10^{4} cm^{-1})')
ylabel('\omega_s  (10^{4} cm^{-1})')
end



%%  For use with fast pump data
%fle = matfile('saved_data/PPcd_with_HEOM_fastpump_with_strong_mode.mat');
fle = matfile('saved_data/PPcd_with_HEOM_fastpump_with_mode_.mat');
%fle = matfile('saved_data/PPcd_with_HEOM_fastpump_no_mode_.mat');

bpset = fle.beam_param_set; ps = 1; t2L = length(t2_range_fs);
cmp = size(bpset{1},3)+ps; tmp = bpset{2}(:,:,ps);
%t2_range_fs = fle.t2_range_fs;
om3_rng = fle.om3_rng;

GSB = imag(fle.R3(:,:,:,cmp)+fle.R4(:,:,:,cmp))/scale_fct2; 
SE = imag(fle.R2(:,:,:,cmp)+fle.R1(:,:,:,cmp))/scale_fct2; 
ESA = imag(fle.R5(:,:,:,cmp)+fle.R6(:,:,:,cmp))/scale_fct2;
GSB2 = imag(fle.R3(:,:,:,cmp+1)+fle.R4(:,:,:,cmp+1))/scale_fct2; 
SE2 = imag(fle.R2(:,:,:,cmp+1)+fle.R1(:,:,:,cmp+1))/scale_fct2; 
ESA2 = imag(fle.R5(:,:,:,cmp+1)+fle.R6(:,:,:,cmp+1))/scale_fct2;

scale_fct2 = max(max(max(abs(fle.R3(:,1,:,1)+fle.R4(:,1,:,1)))));
%% 

%scle= repmat(1./approx_sig_decay,[length(om3_rng),1]);
figure
%sc_fct = 80*10^4;

to_plot = squeeze(0*(SE+ESA+GSB)-(-1)*(SE2+ESA2+GSB2));  %to_plot = to_plot.*sc_fct;
CMRmap2 = interp1(CMRmap,1:0.5:17);
lg = t2_range_fs >= 0 & t2_range_fs< 2000;
contourf(t2_range_fs(lg),om3_rng/1e4,to_plot(:,lg),size(CMRmap,1))
%caxis([-max(abs(to_plot(:))),max(abs(to_plot(:)))])

 colormap(CMRmap2)
 colorbar
xlabel('\tau  (fs)','FontSize',14)
ylabel('\omega_s  (10^{4} cm^{-1})','FontSize',14)

set(gcf, 'Color', 'w');
%%
fig_name = strcat('xxyx_yyz_with_mode.png');
export_fig(fig_name,'-p0.01','-painters')
%  figure
% contourf(t2_range_fs,om3_rng/1e4,squeeze(GSB))
% colorbar
% colormap(CMRmap)
%%
figure
plot(t2_range_fs,to_plot)

%% Compute prony analysis at points

 N_order = 15; lg = t2_range_fs < 1000 & t2_range_fs >= 0; %choose range
t_red= t2_range_fs(lg ); dt2 = (t_red(2)- t_red(1)); fs = 1/dt2;

 Temperature_in_Kelvin = 77;   
[t_scale, B, ~]= inv_cm_unit_sys(Temperature_in_Kelvin);

       pnt = 1.2492*10^4; %1.2255]*10^4;
       %pnt = 12500
    imp_resp = interpn(om3_rng,t_red,ESA(:,lg,1)+SE(:,lg,1)-ESA2(:,lg,1)-SE2(:,lg,1),pnt,t_red);

[num,den] = prony(imp_resp,N_order,N_order);
[r,p,k]= residuez(num,den); %find residues

a_list= r(:); %amplitude list
spoles= log(p(:))*fs; %poles
tau_list= 1./imag(spoles(:));%decay factor
omega_list = imag(spoles(:));
%signals are imag so will be made of cosines
%now collect into sin and cosine factors with imag prefactors
om_new = []; cnt=0; tau_new =[]; sinfct = []; cosfct = [];
[~,neworder]= sort(abs(a_list).*(-tau_list),'descend'); 
a_list = a_list(neworder);
 tau_list = tau_list(neworder);
omega_list = omega_list(neworder);

for lp =1:N_order
    if omega_list(lp) > eps(max(omega_list))
        lpbrk = false;
    for lp2 = [1:lp-1,lp+1:N_order]
        if abs(omega_list(lp2)+omega_list(lp)) < eps(max(omega_list))
           cnt = cnt+1; om_new(cnt) =omega_list(lp); tau_new(cnt) = tau_list(lp);
           cosfct(cnt)= imag(a_list(lp)+a_list(lp2)); 
           sinfct(cnt)= -imag(a_list(lp)-a_list(lp2));
           lpbrk = true;
           break   
        end
    end
    if ~lpbrk %not repeated
       warning('non zero frequency not repeated, assumed to be zero')
       omega_list(lp)
       cnt = cnt+1; om_new(cnt) =0; tau_new(cnt) = tau_list(lp);
           cosfct(cnt)= a_list(lp)/2; 
           sinfct(cnt)= 1i*a_list(lp)/2; 
    end    
        
    elseif abs(omega_list(lp))<= eps(max(omega_list)) %just pure exp decay
        cnt = cnt+1; om_new(cnt) =0;  sinfct(cnt) = 0;
        cosfct(cnt)= imag(a_list(lp)); tau_new(cnt) = tau_list(lp);
    end
end


%{
% Plot the curves
%figure1 = figure('Renderer','painters');
%subplot1 = subplot(2,1,1,'Parent',figure1);
%subplot2 = subplot(2,1,2,'Parent',figure1);
%}
 curve = zeros(cnt,length(t_red));
    for lp = 1:cnt
        curve(lp,:) = cosfct(lp).*cos(om_new(lp)*(t_red-t_red(1))) +...
                      sinfct(lp).*sin(om_new(lp)*(t_red-t_red(1)))   ;
        curve(lp,:) = curve(lp,:).*exp(t_red./tau_new(lp)) ;        
    end
    %{
%     subplot(2,1,1);
%     plot1 = plot(t_red,imag(curve),'Parent',subplot1);
% 
%     xlabel('t_2 (fs)','FontSize',16); ylabel('ES signal, a.u.','FontSize',14);
%     for lp =1:cnt
%         set(plot1(lp),'DisplayName',strcat('\omega=',num2str(om_new(lp)*1e3/t_scale,3),'cm^{-1}'))
%     end
% legend('show')
% 
% plot(t_red,[sum(curve);imp_resp],'Parent',subplot2);
    %}
prony_compare_fig(t_red,imag(curve(1:4,:)), [sum(curve);imp_resp],om_new*1e3/t_scale)
%%
set(gcf, 'Color', 'w');
fig_name = strcat('prony_compare_12492_yyyx_with_mode.png');
export_fig(fig_name,'-p0.01','-painters')