%fle = matfile('saved_data/plenio_cd_config_with_RF_nomode.mat');
fle = matfile('saved_data/plenio_cd_config_with_RF3.mat');

%%
om1_rng = [12320,12520]; 
om3_rng = linspace(1.21e4,1.27e4,300);

scale_fct2 = max(max(max(abs(fle.R3(:,:,:,1)+fle.R4(:,:,:,1)))));
to_plot1 = zeros(length(om3_rng),length(t2_range_fs),2);
to_plot2 = to_plot1;  to_plot3 = to_plot1; 
cmp = 3;
for lp=0:1
to_plot1(:,:,lp+1) = imag(fle.R3(:,:,length(om1_rng)+1,cmp+lp)+fle.R4(:,:,length(om1_rng)+1,cmp+lp))/scale_fct2; 
to_plot2(:,:,lp+1)  = imag(fle.R1(:,:,length(om1_rng)+1,cmp+lp)+fle.R2(:,:,length(om1_rng)+1,cmp+lp))/scale_fct2; 
to_plot3(:,:,lp+1)  = imag(fle.R5(:,:,length(om1_rng)+1,cmp+lp)+fle.R6(:,:,length(om1_rng)+1,cmp+lp))/scale_fct2;
end   
%%
figure
pcolor(t2_range_fs,om3_rng,abs(to_plot2(:,:,1)+to_plot3(:,:,1)...
                            -to_plot2(:,:,2)-to_plot3(:,:,2)))
shading flat
figure
pcolor(t2_range_fs,om3_rng,abs(to_plot2(:,:,2)+to_plot3(:,:,2)))
shading flat
%%
figure
pcolor(t2_range_fs,om3_rng,(to_plot2(:,:,1)+to_plot3(:,:,1)...
                            -to_plot2(:,:,2)-to_plot3(:,:,2)))
shading flat
figure
pcolor(t2_range_fs,om3_rng,(to_plot2(:,:,2)+to_plot3(:,:,2)))
shading flat


%%  prony decomp

%om_prony = [12320,12520];
%om_prony = [12310,12320,12330,12510,12520,12530];
om_prony = [12260,12320,12510,12550];
t_red = t2_range_fs(t2_range_fs<500)/1e3*t_scale; 
 dt2 = (t_red(2)- t_red(1)); fs = 1/dt2;

N_order = 9;
a_list = complex(zeros(length(om_prony),N_order)); spoles = a_list;
tau_list = zeros(length(om_prony),N_order); omega_list = tau_list;
% figure
clear om_pos a_sin a_cos

for lpval = 1:length(om_prony)

    [~,tmp] = min(abs(om3_rng-om_prony(lpval)));
    
impulse_resp = (to_plot2(tmp,:,1)+to_plot3(tmp,:,1)...
                -0*(to_plot2(tmp,:,2)+to_plot3(tmp,:,2))).';
[num,den] = prony(impulse_resp,N_order,N_order);
[r,p,k]= residuez(num,den); %find residues

a_list(lpval,:)= r(:); %amplitude list
spoles(lpval,:)= log(p(:))*fs; %poles
tau_list(lpval,:)= 1./real(spoles(lpval,:));%decay factor
omega_list(lpval,:) = imag(spoles(lpval,:));

om_pos{lpval} = omega_list(lpval,omega_list(lpval,:)>=0);
        a_sin{lpval} = imag(a_list(lpval,omega_list(lpval,:)>=0));
        a_cos{lpval} = real(a_list(lpval,omega_list(lpval,:)>=0));
end

%% Assume made of sines and cosines and take freq out

figure
stem(om_pos.',log(a_sin.'))
figure
stem(om_pos.',log(a_cos.'))


%% Isolate low frequency contributions
colour_spread = jet;
colour_spread = colour_spread(1:length(colour_spread)/length(om_prony)...
                                :length(colour_spread),:);

figure
hold on
for lpval = 1:length(om_prony)
lg = abs(omega_list(lpval,:)) < 100;
tmp = find(lg);
curve = zeros(length(tmp),length(t_red));
for lp2 = 1:length(tmp)
curve(lp2,:) = a_list(lpval,tmp(lp2)).*exp(spoles(lpval,tmp(lp2))*t_red);
end
plot(t2_range_fs(t2_range_fs<500),sum(curve,1),'color',colour_spread(lpval,:))
end

%% For frequency resolution
%fle = matfile('saved_data/PPcd_with_HEOM_with_mode.mat');
fle = matfile('saved_data/PPcd_with_HEOM_no_mode.mat');

 om_pnts=[12320,12320;12320,12520;12520,12320;12520,12520];
om1_rng = linspace(1.20e4,1.28e4,50); %range in cm^(-1)
om3_rng = linspace(1.20e4,1.28e4,150) ;
[a,b] = min(abs(om1_rng-12320)); [a,b2] = min(abs(om1_rng-12520));
[a,c] = min(abs(om3_rng-12320)); [a,c2] = min(abs(om3_rng-12520));

%%
cmp = 3;
t2_range_fs = linspace(0,2000,210); 
scale_fct2 = max(max(max(abs(fle.R3(:,:,:,1)+fle.R4(:,:,:,1)))));
to_plot1 = zeros(length(t2_range_fs),size(om_pnts,1),2);
to_plot2 = to_plot1;  to_plot3 = to_plot1; 

        to_plot1(:,1,:) = fle.R3(c,:,b,cmp:cmp+1)+fle.R4(c,:,b,cmp:cmp+1);
        to_plot1(:,2,:) = fle.R3(c2,:,b,cmp:cmp+1)+fle.R4(c2,:,b,cmp:cmp+1);       
        to_plot1(:,3,:) = fle.R3(c2,:,b2,cmp:cmp+1)+fle.R4(c2,:,b2,cmp:cmp+1);
        to_plot1(:,4,:) = fle.R3(c,:,b2,cmp:cmp+1)+fle.R4(c,:,b2,cmp:cmp+1); 
   
        to_plot2(:,1,:) = fle.R1(c,:,b,cmp:cmp+1)+fle.R2(c,:,b,cmp:cmp+1);
        to_plot2(:,2,:) = fle.R1(c2,:,b,cmp:cmp+1)+fle.R2(c2,:,b,cmp:cmp+1);       
        to_plot2(:,3,:) = fle.R1(c2,:,b2,cmp:cmp+1)+fle.R2(c2,:,b2,cmp:cmp+1);
        to_plot2(:,4,:) = fle.R1(c,:,b2,cmp:cmp+1)+fle.R2(c,:,b2,cmp:cmp+1); 
   
        to_plot3(:,1,:) = fle.R5(c,:,b,cmp:cmp+1)+fle.R6(c,:,b,cmp:cmp+1);
        to_plot3(:,2,:) = fle.R5(c2,:,b,cmp:cmp+1)+fle.R6(c2,:,b,cmp:cmp+1);       
        to_plot3(:,3,:) = fle.R5(c2,:,b2,cmp:cmp+1)+fle.R6(c2,:,b2,cmp:cmp+1);
        to_plot3(:,4,:) = fle.R5(c,:,b2,cmp:cmp+1)+fle.R6(c,:,b2,cmp:cmp+1); 

to_plot1 = to_plot1/scale_fct2;  to_plot2 = to_plot2/scale_fct2; 
 to_plot3 = to_plot3/scale_fct2; 

 %%
figure
plot0 = plot(t2_range_fs,imag(to_plot1(:,:,1)));
xlabel('t_2 (fs)'); ylabel('GS signal, a.u.');
for lp=1:4
    set(plot0(lp),'DisplayName',strcat('[om_1,om_3]=[',num2str(om_pnts(lp,1)),...
                        ',',num2str(om_pnts(lp,2)),']cm^{-1}'));
end
set(plot0(lp),'Color',[0.749019622802734 0 0.749019622802734]);
legend('show');
%%
figure
plot1 = plot(t2_range_fs,to_plot2(:,:,1)+to_plot3(:,:,1)-0*(...
    to_plot2(:,:,2)+to_plot3(:,:,2)));
xlabel('t_2 (fs)'); ylabel('ES signal, a.u.');
for lp=1:4
    set(plot1(lp),'DisplayName',strcat('[om_1,om_3]=[',num2str(om_pnts(lp,1)),...
                        ',',num2str(om_pnts(lp,2)),']cm^{-1}'));
end
set(plot1(lp),'Color',[0.749019622802734 0 0.749019622802734]);
legend('show');
 %%  prony decomp

%om_prony = [12320,12520];
%om_prony = [12310,12320,12330,12510,12520,12530];
om_prony = [12260,12320,12510,12550];
t_red = t2_range_fs(t2_range_fs<1000)/1e3*t_scale; 
 dt2 = (t_red(2)- t_red(1)); fs = 1/dt2;

N_order = 9;
a_list = complex(zeros(length(om_prony),N_order)); spoles = a_list;
tau_list = zeros(length(om_prony),N_order); omega_list = tau_list;
% figure
clear om_pos a_sin a_cos

for lpval = 1:4
impulse_resp = imag(to_plot1(:,lpval,1)) ;   
%impulse_resp = imag((to_plot2(:,lpval,1)+to_plot3(:,lpval,1)...
%                -0*(to_plot2(:,lpval,2)+to_plot3(:,lpval,2))).');
[num,den] = prony(impulse_resp,N_order,N_order);
[r,p,k]= residuez(num,den); %find residues

a_list(lpval,:)= r(:); %amplitude list
spoles(lpval,:)= log(p(:))*fs; %poles
tau_list(lpval,:)= 1./real(spoles(lpval,:));%decay factor
omega_list(lpval,:) = imag(spoles(lpval,:));

om_pos{lpval} = omega_list(lpval,omega_list(lpval,:)>=-10);
        a_sin{lpval} = imag(a_list(lpval,omega_list(lpval,:)>=-10));
        a_cos{lpval} = real(a_list(lpval,omega_list(lpval,:)>=-10));
end




%% Isolate low frequency contributions
% colour_spread = jet;
% colour_spread = colour_spread(1:length(colour_spread)/length(om_prony)...
%                                 :length(colour_spread),:);
colour_spread = [1,0,0;0,1,0;0,0,1;0.749019622802734 0 0.749019622802734];
figure1 = figure('Renderer','painters','InvertHardcopy','off',...
    'Color',[1 1 1]);
axes1 = axes('Parent',figure1,'FontSize',14);
hold on
for lpval = 1:length(om_prony)
lg = abs(omega_list(lpval,:)) < 100;
tmp = find(lg);
curve = zeros(length(tmp),length(t_red));
for lp2 = 1:length(tmp)
curve(lp2,:) = a_list(lpval,tmp(lp2)).*exp(spoles(lpval,tmp(lp2))*t_red);
end
plot(t2_range_fs(t2_range_fs<1000),real(sum(curve,1)),'color',colour_spread(lpval,:),...
        'DisplayName',strcat('[om_1,om_3]=[',num2str(om_pnts(lpval,:)),']cm^{-1}'),...
            'Parent',axes1,'LineWidth',2)

end
xlabel('t_2 (fs)','FontSize',16); ylabel('ES signal, a.u.','FontSize',14);