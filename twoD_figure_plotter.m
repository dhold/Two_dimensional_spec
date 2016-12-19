%% Purely linear spectra
scale_fct = max(real(lin_spec{1}{1})); 
 figure
 plot(om1_rng,real(lin_spec{1}{1})/scale_fct);
 hold on
 plot(om3_rng,real(lin_spec{2}{1})/scale_fct);
%%
figure
plotyy(om1_rng,real(lin_spec{1}{1})/scale_fct,om1_rng,imag(lin_spec{1}{2})/scale_fct)
hold on
tmp1 = interp1(om3_rng,real(lin_spec{2}{1}),om1_rng,'pchip');
tmp2 = interp1(om3_rng,imag(lin_spec{2}{2}),om1_rng,'pchip');
plotyy(om1_rng,tmp1/scale_fct,om1_rng,tmp2/scale_fct)
%%
figure
plotyy(om1_rng,imag(lin_spec{1}{1})/scale_fct,om1_rng,real(lin_spec{1}{2})/scale_fct)
tmp1 = interp1(om3_rng,imag(lin_spec{2}{1}),om1_rng,'pchip');
tmp2 = interp1(om3_rng,real(lin_spec{2}{2}),om1_rng,'pchip');
hold on
plotyy(om1_rng,tmp1/scale_fct,om1_rng,tmp2/scale_fct)



%% PP (only really valid when k1=k2 = pol1=pol2)
cmp =4;
tmp = squeeze(real(R3(:,1,:,cmp)+R4(:,1,:,cmp)));
%upper_selec = om3_rng > 1.24e4;
[a,b] = max(abs(tmp),[],cmp); [a2,b2] =  max(tmp(b(1),:));% b2 = 20
tmp2 = real(squeeze(R2(:,:,b2,cmp)+R1(:,:,b2,cmp))); %SE
%tmp2 =  real(R3(:,:,b2,cmp)+R4(:,:,b2,cmp)); %GSB
%tmp2 =  real(squeeze(R5(:,:,b2,cmp)+R6(:,:,b2,cmp))); %ESA
figure
pcolor(t2_range,om3_rng,tmp2)
shading flat
figure
plot(t2_range,tmp2(max(b(1)-10,1):2:b(1)+5,:))
%% Photon echo at point in time
pnt = length(t2_range); 
cmp = 3;

figure
pcolor(om3_rng,om1_rng,(imag(squeeze(R2(:,pnt,:,cmp)+R1(:,pnt,:,cmp)))).')
xlabel('\omega_3')
ylabel('\omega_1')
shading flat
figure
pcolor(om3_rng,om1_rng,(imag(squeeze(R3(:,pnt,:,cmp)+R4(:,pnt,:,cmp))).'))
shading flat
xlabel('\omega_3')
ylabel('\omega_1')
figure
pcolor(om3_rng,om1_rng,(imag(squeeze(R5(:,pnt,:,cmp)+ R6(:,pnt,:,cmp))).'))
shading flat
xlabel('\omega_3')
ylabel('\omega_1')

%% Pump probe signal
comp = 3; tau_fs = 1;  %pulse width
%tau_cm = tau_fs * t_scale / 1000; dw = 1/tau_cm ;
tau_cm = 1; 

om0_probe = 1001;

E_probe = @(om) exp(-(om-om0_probe).^2*tau_cm^2/2)*tau_cm; %no need to normalise properly

PP_sig = squeeze(R2(:,:,:,comp)+R1(:,:,:,comp)+R3(:,:,:,comp)...
                 +R4(:,:,:,comp)+R5(:,:,:,comp)+R6(:,:,:,comp));
PP_sig = repmat(reshape(om3_rng.*E_probe(om3_rng),[length(om3_rng),1])...
        ,[1,length(t2_range)]) .*trapz(om1_rng,PP_sig,3)/2/pi;

figure
pcolor(t2_range,om3_rng,real(PP_sig))
shading flat
xlabel('\tau ')
ylabel('\omega_1')
    
    
    
%% total rephasing and non rephasing
pnt = 2; comp=1;

figure
pcolor(om1_rng,om3_rng,(real(squeeze(R2(:,pnt,:,comp) + R3(:,pnt,:,comp)+R5(:,pnt,:,comp)))))
shading interp
xlabel('\omega_1')
ylabel('\omega_3')

figure
pcolor(om1_rng,om3_rng,(real(squeeze(R1(:,pnt,:,comp) + R4(:,pnt,:,comp)+R6(:,pnt,:,comp)))))
shading interp
xlabel('\omega_1')
ylabel('\omega_3')
%% nonRephasing
pnt=1; cmp = 4;
figure
pcolor(om1_rng,om3_rng,(real(squeeze(R1(:,pnt,:,cmp)))))
shading interp
xlabel('\omega_1')
ylabel('\omega_3')

figure
pcolor(om1_rng,om3_rng,(real(squeeze(R4(:,pnt,:,cmp)))))
shading interp
xlabel('\omega_1')
ylabel('\omega_3')

figure
pcolor(om1_rng,om3_rng,(real(squeeze(R6(:,pnt,:,cmp)))))
shading interp
xlabel('\omega_1')
ylabel('\omega_3')
%% Rephasing

scale_fct2 = max(max(max(abs(real(R3(:,:,:,1))))));
pnt = 1%length(t2_range);
figure
pcolor(om3_rng,om1_rng,abs(real(squeeze(R3(:,pnt,:,1)))).'/scale_fct2)
shading interp
xlabel('\omega_3')
ylabel('\omega_1')

%%
figure
pcolor(om3_rng,om1_rng,abs(real(squeeze(R3(:,pnt,:,1)+R2(:,pnt,:,1)+R5(:,pnt,:,1))).')/scale_fct2)
shading interp
xlabel('\omega_3')
ylabel('\omega_1')

%%
figure
pcolor(om1_rng,om3_rng,abs(real(squeeze(R2(:,pnt,:,1))))/scale_fct2)
shading interp
%%
figure
pcolor(om1_rng,om3_rng,abs(real(squeeze(R5(:,pnt,:,1))))/scale_fct2)
shading interp

%% This one is used if it is in time domain
%{
tmp4 = fftshift(ifft(fftshift(ifft(squeeze(R4(:,7,:,1)),[],1),1),[],2),2);
tmp1 = fftshift(ifft(fftshift(ifft(squeeze(R1(:,7,:,1)),[],1),1),[],2),2);
tmp6 = fftshift(ifft(fftshift(ifft(squeeze(R6(:,7,:,1)),[],1),1),[],2),2);
figure
pcolor(om1 ,om3 ,real(tmp1))
shading flat
figure
pcolor(om1 ,om3 ,real(tmp4))
shading flat
figure
pcolor(om1 ,om3 ,real(tmp6))
shading flat
%}