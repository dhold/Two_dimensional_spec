syms a b c t1 t2 t3 t4
ry = [cos(b),0,sin(b);0,1,0;-sin(b),0,cos(b)];
rz = [cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
rz2 = [cos(c),-sin(c),0;sin(c),cos(c),0;0,0,1];

rot_op = rz2*ry*rz;

e3 = subs(rot_op*[1;0;0],[c,a],[pi/4,t3]);
e4 = subs(rot_op*[1;0;0],[c,a],[3*pi/4,t4]);
e1 = subs(rot_op*[1;0;0],[c,a],[7*pi/4,t1]);
e2 = subs(rot_op*[1;0;0],[c,a],[5*pi/4,t2]);

k1 = subs(rot_op*[0;0;1],[c,a],[7*pi/4,t1]); %+kx - ky
k2 = subs(rot_op*[0;0;1],[c,a],[5*pi/4,t2]); %-kx-ky
k3 = subs(rot_op*[0;0;1],[c,a],[pi/4,t3]); %+kx+ky
k4 = subs(rot_op*[0;0;1],[c,a],[3*pi/4,t4]); %-kx + ky

test = [e1.'*e2 * e4.'*e3,e1.'*e3 * e2.'*e4, e1.'*e4 * e2.'*e3];
test= test*[4,-1,-1;-1,4,-1;-1,-1,4]/30; test = test.';

%cannot find general analytic solution
%test2 = solve(test == 0,t1,t2,t3);
%test3 = solve([test;t4] == [0;0;0;-3*pi/4],t1,t2,t3); % not over determined
%{
%tset = [0,0,0,0 ; pi/4*[1,3,7,5];-pi/4*[1,3,7,5]];
%tset = pi/4*[2,2,-2,-2 ; 1,0,1,0 ];

%    e1.'*e4 .* e2.' *e3==0,t1,t2,t3,t4);
%
% tmp = subs([e1.'*e2 .* e3.' *e4,e1.'*e3 .* e2.' *e4,e1.'*e4 .* e2.' *e3],b,0);
% tmp = simplify(tmp);
% soln = solve(tmp==[0,0,0],t1,t2,t3,t4);
% tmp = taylor([e1.'*e2 .* e3.' *e4,e1.'*e3 .* e2.' *e4,e1.'*e4 .* e2.' *e3],b,'Order',1); 
% soln2 = solve(tmp==[0,0,0],t1,t2,t3,t4);

%[x,fval,exitflag] = fminunc(ff ,rand([4,1]));
% fval_save = zeros(10000,1); x_save =  zeros(10000,4);
% for lp = 1:10000
%     tmp = 2*pi*rand(4,1); x_save (lp,:) = tmp;
%     fval_save (lp) = ff(tmp);
% end
%}

test2 = solve(e1.'*e2 ==0, e1.'*e3==0,e1.'*e4==0,t1,t2,t3);
%% Minimisation stuff

b_range = pi*linspace(0,1/32,100);

OPTIONS.InitTemp = 1e-8; OPTIONS.CoolSched = @(T) 0.9.*T;
OPTIONS.StopTemp = 1e-10; OPTIONS.StopVal = 1e-8;
OPTIONS.MaxConsRej = 4000;  OPTIONS.MaxTries = 1000;
%{
% for lp = 1:length(b_range)
%     
%     tmp = subs(abs(e1.'*e2 .* e3.' *e4) + abs(e1.'*e3 .* e2.' *e4)+...
%     abs(e1.'*e4 .* e2.' *e3),'b',b_range(lp));
% f = matlabFunction(tmp);
% ff = @(x) f(x(1),x(2),x(3),x(4));
%     if lp > 1; guess = minimum{lp-1};
%     else; guess = pi/4*[-1,-1,3,-3]; end %pi/4*[1,-1,3,3];
% [minimum{lp},fval{lp}] = anneal(ff,guess,OPTIONS);
% OPTIONS.InitTemp = OPTIONS.InitTemp*2;
% end
%}
%{
clear minimum fval minimum2 fval2  minimum3 fval3 
options.TolFun = 1e-8;  options.TolX = 1e-8; 
cnst = 3.004296875000003;
for lp = 1:length(b_range)
    
    tmp = subs(abs(e1.'*e2 .* e3.' *e4) + abs(e1.'*e3 .* e2.' *e4)+...
    abs(e1.'*e4 .* e2.' *e3),'b',-b_range(lp));
f = matlabFunction(tmp);
init_rot = pi/4*[-1,-1,3,-3];
ff = @(x) f(x(1)+init_rot(1),x(2)+init_rot(2),x(3)+init_rot(3),x(4)+init_rot(4));
ff2 = @(x) f(x(1)+init_rot(1),x(2)+init_rot(2),-x(1)+init_rot(3),-x(2)+init_rot(4));
ff3 = @(x) f(x+init_rot(1),x*cnst+init_rot(2),-x+init_rot(3),-x*cnst+init_rot(4));
    if lp > 1; guess3 = minimum3(lp-1); %guess = minimum{lp-1}; guess2 = minimum2{lp-1}; 
    else; guess = [0,0,0,0]; guess2 = [0,0]; guess3=0; end
%[minimum{lp},fval{lp}] =  fminunc(ff ,guess);    %more general
%[minimum2{lp},fval2{lp}] =  fminunc(ff2 ,guess2);  %less general
[minimum3(lp),fval3(lp)] =  fminunc(ff3 ,guess3,options);  %single param fir
end
%}
clear minimum fval minimum2 fval2  minimum3 fval3 
options.TolFun = 1e-8;  options.TolX = 1e-8; 
cnst = 3;
for lp = 1:length(b_range)
    
    %tmp = subs([(e1.'*e2 .* e3.' *e4) , (e1.'*e3 .* e2.' *e4),...
     %   (e1.'*e4 .* e2.' *e3)],'b',-b_range(lp));
     tmp = subs([(e1.'*e2 .* e3.' *e4) , (e1.'*e3 .* e2.' *e4),...
        (e1.'*e4 .* e2.' *e3)],'b',b_range(lp)); %v slow step actually
    tmp2 = 1/30*tmp*[4,-1,-1;-1,4,-1;-1,-1,4] ; tmp3 = sum(abs(tmp2));
    
f = matlabFunction(tmp3);
%init_rot = pi/4*[3,3,-1,-3]; %chiral pump
init_rot = pi/4*[1,1,-1,-3]; %chiral pump2
%init_rot =pi/4*[1,3,1,-3]; %chiral probe
ff = @(x) f(x(1)+init_rot(1),x(2)+init_rot(2),x(3)+init_rot(3),x(4)+init_rot(4));
ff2 = @(x) f(x(1)+init_rot(1),-x(1)+init_rot(2),x(2)+init_rot(3),-x(2)+init_rot(4));
ff3 = @(x) f(x+init_rot(1),-x+init_rot(2),x*cnst+init_rot(3),-x*cnst+init_rot(4));
    if lp > 1; guess3 = minimum3(lp-1); guess = minimum{lp-1}; guess2 = minimum2{lp-1}; 
    else; guess = [0,0,0,0]; guess2 = [0,0]; guess3=0; end
[minimum{lp},fval{lp}] =  fminunc(ff ,guess);    %more general
[minimum2{lp},fval2{lp}] =  fminunc(ff2 ,guess2);  %less general
[minimum3(lp),fval3(lp)] =  fminunc(ff3 ,guess3,options);  %single param fir
[testtt(lp,:) ,testtt2(lp)]= fminunc(ff2 ,[1,3]/4);
end

%% Analytic estimates chiral probe

init_rot =pi/4*[1,3,1,-3];  %chiral probe
b_range = pi*linspace(0,1/6,1000);
ff_tmp = matlabFunction([(e1.'*e2 .* e3.' *e4),(e1.'*e3 .* e2.' *e4),...
    e1.'*e4 .* e2.' *e3]);

ff_tmp2 = @(b) ff_tmp(b,-(3*b^2/4)+init_rot(1),(3*b^2/4)+init_rot(2),...
                    -b^2/4+init_rot(3),b^2/4+init_rot(4));
%ff_tmp3 = @(b) ff_tmp(b,-(3*b^2/4+0.62744*b^4)+init_rot(1),(3*b^2/4+0.62744*b^4)+init_rot(2),...
%           -(b^2/4+0.045117*b^4)+init_rot(3),(b^2/4+0.045117*b^4)+init_rot(4)); 
ff_tmp3 = @(b) ff_tmp(b,-(3*b^2/4+15*b^4/24)+init_rot(1),(3*b^2/4+15*b^4/24)+init_rot(2),...
           -(b^2/4+b^4/24)+init_rot(3),(b^2/4+b^4/24)+init_rot(4)); 

 fb_l = @(b) b.^2/4+b.^4/24 + (0.6840*b.^6  -1.7634*b.^8-0.6270*b.^10)*1e-3; 
 fb_u = @(b) 3*b.^2/4+15*b.^4/24 + 0.4146*b.^6  +  0.2183*b.^8  -0.3670*b.^10;       
       ff_tmp_test = @(b) ff_tmp(b,-fb_u(b)+init_rot(1),fb_u(b)+init_rot(2),...
           - fb_l(b)+init_rot(3), fb_l(b)+init_rot(4));           
 for k = 1:length(b_range); 
     toplot0(k) = sum(abs(ff_tmp(b_range(k),init_rot(1),init_rot(2),init_rot(3),init_rot(4))));
     toplot1(k) = sum(abs(ff_tmp2(b_range(k))));
    toplot2(k) = sum(abs(ff_tmp3(b_range(k))));
    toplot3(k) = sum(abs(ff_tmp_test(b_range(k))));
 end

%rnj = cat(1,fval2{:})< 3e-7;        
%       tmp = reshape([minimum2{:}],2,[]);
%tmp2 = ezfit( b_range(rnj),tmp(1,rnj),'-3*x^2/4-15*x^4/24+a*x^6');
%tmp3 = ezfit( b_range(rnj),tmp(2,rnj),'-x^2/4-x^4/24+a*x^6');

figure; plot(b_range, [toplot0;toplot1;toplot2;toplot3])

%% Plot angles between beams
b_range = pi*linspace(0,1/6,1000);
 fb_l = @(b) b.^2/4+b.^4/24 + (0.6840*b.^6  -1.7634*b.^8-0.6270*b.^10)*1e-3; 
 fb_u = @(b) 3*b.^2/4+15*b.^4/24 + 0.4146*b.^6  +  0.2183*b.^8  -0.3670*b.^10;    

 figure1 = figure;
axes1 = axes('Parent',figure1);  hold(axes1,'on');

inter_angle = zeros(length(b_range),6);
 for k = 1:length(b_range);
     b = b_range(k);
     inter_angle(k,1) = acos(dot(-fb_u(b)+init_rot(1),fb_u(b)+init_rot(2))); %12
     inter_angle(k,2) = acos(dot(-fb_u(b)+init_rot(1),-fb_l(b)+init_rot(3))); %13
     inter_angle(k,3) = acos(dot(-fb_u(b)+init_rot(1),fb_l(b)+init_rot(4))); %14
     inter_angle(k,4) = acos(dot(fb_u(b)+init_rot(2),-fb_l(b)+init_rot(3))); %23
     inter_angle(k,5) = acos(dot(fb_u(b)+init_rot(2),fb_l(b)+init_rot(4))); %24
     inter_angle(k,6) = acos(dot(-fb_l(b)+init_rot(3),fb_l(b)+init_rot(4))); %34
 end      
plot1 = plot(b_range,inter_angle-repmat(inter_angle(1,:),length(b_range),1),'Parent',axes1);
set(plot1(1),'DisplayName','\theta_{12}');
set(plot1(2),'DisplayName','\theta_{13}');
set(plot1(3),'DisplayName','\theta_{14}','LineStyle','--');
set(plot1(4),'DisplayName','\theta_{23}');
set(plot1(5),'DisplayName','\theta_{24}','LineStyle',':');
set(plot1(6),'DisplayName','\theta_{34}','LineStyle','-.');

xlabel('\Theta (radians)');  % xlim(axes1,[0 0.5]);
ylabel('Angle between polarization (radians)');

box(axes1,'on');
legend1 = legend(axes1,'show');    set(legend1,'Position',...
[0.138834329030446 0.485804414711539 0.17092337735978 0.424815971329887]);      
       
%% Analytic estimates chiral pump

options.TolFun = 1e-8;  options.TolX = 1e-8; 
init_rot = pi/4*[3,3,-1,-3]; %chiral pump
%init_rot = pi/4*[1,1,-1,-3]; %other chiral pump
ff_tmp = matlabFunction([(e1.'*e2 .* e3.' *e4),(e1.'*e3 .* e2.' *e4),...
    e1.'*e4 .* e2.' *e3]);

ff_tmp2 = @(b) ff_tmp(b,(b^2/4)+init_rot(1),-(b^2/4)+init_rot(2),...
           (3*b^2/4)+init_rot(3),-(3*b^2/4)+init_rot(4));
ff_tmp3 = @(b) ff_tmp(b,(b^2/4+b^4/24)+init_rot(1),-(b^2/4+b^4/24)+init_rot(2),...
           (3*b^2/4+15*b^4/24)+init_rot(3),-(3*b^2/4+15*b^4/24)+init_rot(4)); 
ff_tmp4 = @(c,b) ff_tmp(b,(b^2/4+b^4/24+c(1))+init_rot(1),-(b^2/4+b^4/24+c(1))+init_rot(2),...
           (3*b^2/4+15*b^4/24+c(2))+init_rot(3),-(3*b^2/4+15*b^4/24+c(2))+init_rot(4));     
 
 fb_l = @(b) b.^2/4+b.^4/24 + (0.6840*b.^6  -1.7634*b.^8-0.6270*b.^10)*1e-3; 
 fb_u = @(b) 3*b.^2/4+15*b.^4/24 + 0.4146*b.^6  +  0.2183*b.^8  -0.3670*b.^10;       
ff_tmp_test = @(b) ff_tmp(b, fb_l(b)+init_rot(1),-fb_l(b)+init_rot(2),...
            fb_u(b)+init_rot(3),- fb_u(b)+init_rot(4)); 
       %T1 = -T2; T3 = -T4;
%T1 = t2^2/2 + t2^4/24; T3 = 3*t2^2/2 + 15*t2^4/24;
b_range = pi*linspace(0,1/6,1000);
tmp = zeros(3,1000); tmp2 = tmp; tmp4 = tmp;
tmp3 = zeros(1,1000); cmin = zeros(2,1000);
for lp = 1:1000
    tmp(:,lp) =  ff_tmp2(b_range(lp)); tmp2(:,lp) =  ff_tmp3(b_range(lp));
    tmp4(:,lp) = ff_tmp_test(b_range(lp));
    ff_tmp5 = @(c) sum(abs(ff_tmp4(c,b_range(lp))));
    if lp > 1; guess = cmin(:,lp-1); 
    else; guess = [0;0]; end    
   % [cmin(:,lp),tmp3(lp)] =  fminunc(ff_tmp5 ,guess);
end
%HOT1 = ezfit(b_range(tmp3<1e-7),cmin(1,tmp3<1e-7),'a*x^6+b*x^8+c*x^10');
%HOT2 = ezfit(b_range(tmp3<1e-7),cmin(2,tmp3<1e-7),'a*x^6+b*x^8+c*x^10');
figure; plot(b_range, [sum(abs(tmp));sum(abs(tmp2));sum(abs(tmp4))])
%% Calculate the strength of the chiral bits

%ee3 = subs(rot_op*[1;0;0],[c,a],[pi/4,(3*b^2/4+15*b^4/24+b^6*0.40312)+init_rot(3)]);
%ee4 = subs(rot_op*[1;0;0],[c,a],[3*pi/4,-(3*b^2/4+15*b^4/24+b^6*0.40312)+init_rot(4)]);
ee3 = subs(rot_op*[1;0;0],[c,a],[pi/4,(3*b^2/4+15*b^4/24)+init_rot(3)]);
ee4 = subs(rot_op*[1;0;0],[c,a],[3*pi/4,-(3*b^2/4+15*b^4/24)+init_rot(4)]);
ee1 = subs(rot_op*[1;0;0],[c,a],[7*pi/4,(b^2/4+b^4/24)+init_rot(1)]);
ee2 = subs(rot_op*[1;0;0],[c,a],[5*pi/4,-(b^2/4+b^4/24)+init_rot(2)]);

FF1 = matlabFunction(ori_F_calc([ee1,ee2,ee3,ee4,k1].'));
FF2 = matlabFunction(ori_F_calc([ee1,ee2,ee3,ee4,k2].'));
FF3 = matlabFunction(ori_F_calc([ee1,ee2,ee3,ee4,k3].'));
FF4 = matlabFunction(ori_F_calc([ee1,ee2,ee3,ee4,k4].'));
FF = @(b) sum(FF1(b).^2)+sum(FF2(b).^2)+sum(FF3(b).^2)+sum(FF4(b).^2);
Ff = @(b) FF1(b).^2+FF2(b).^2+FF3(b).^2+FF4(b).^2;

GG1 = matlabFunction(ori_F_calc([cross(ee1,k1),ee2,ee3,ee4].'));
GG2 = matlabFunction(ori_F_calc([ee1,cross(ee2,k2),ee3,ee4].'));
GG3 = matlabFunction(ori_F_calc([ee1,ee2,cross(ee3,k3),ee4].'));
GG4 = matlabFunction(ori_F_calc([ee1,ee2,ee3,cross(ee4,k4)].'));
tmp=[taylor(ori_F_calc([cross(ee1,k1),ee2,ee3,ee4].'),'b'),taylor(ori_F_calc([ee1,cross(ee2,k2),ee3,ee4].')),...
 taylor(ori_F_calc([ee1,ee2,cross(ee3,k3),ee4].'),'b'),taylor(ori_F_calc([ee1,ee2,ee3,-cross(ee4,k4)].'))];   

Gf = @(b) GG1(b).^2+GG2(b).^2+GG3(b).^2+GG4(b).^2;
%%
figure; 
subplot(2,2,1); plot(b_range,GG1(b_range) - repmat(GG1(b_range(1)),[1,length(b_range)]));
subplot(2,2,2); plot(b_range,GG2(b_range) - repmat(GG2(b_range(1)),[1,length(b_range)]));
subplot(2,2,3); plot(b_range,GG3(b_range) - repmat(GG3(b_range(1)),[1,length(b_range)]));
subplot(2,2,4); plot(b_range,GG4(b_range) - repmat(GG4(b_range(1)),[1,length(b_range)]));
%% Single parameter minimisation analysis
ff4= matlabFunction([(e1.'*e2 .* e3.' *e4),(e1.'*e3 .* e2.' *e4),...
    e1.'*e4 .* e2.' *e3]);
%tmp = ezfit( b_range,minimum3,'a*x+b*x^2+c*x^3+d*x^4');
tmp = ezfit( b_range(1:50),minimum3(1:50),'a*x+b*x^2');%+d*x^4');
%angle_fn =  tmp.m(1).*b_range+tmp.m(2).*b_range.^2;
pp = splinefit(b_range,minimum3,3); %interpolating spline fit
angle_fn2 = ppval(pp,b_range);
%{
%angle_fn(1:50) =  tmp.m(1).*b_range(1:50)+tmp.m(2).*b_range(1:50).^2;
%tmp2 = ezfit( b_range(50:end),minimum3(50:end)-angle_fn(50),'a*x+b*x^2');%+d*x^4');
%angle_fn(51:300) =  tmp2.m(1).*b_range(51:end)+tmp2.m(2).*b_range(51:end).^2+angle_fn(50);
          %  tmp.m(3).*b_range.^3+tmp.m(4).*b_range.^4;
%tmp = ezfit( b_range.',minimum3,'a*x^2+b*x');
%angle_fn =  tmp.m(1).*b_range.^2+tmp.m(2).*b_range;
%}
min_smooth = zeros(3,length(b_range)); min_act = min_smooth;
for lp = 1:length(b_range)
    min_smooth(:,lp) = ff4(b_range(lp),init_rot(1)+angle_fn2(lp),...
        init_rot(2)-angle_fn2(lp),init_rot(3)+cnst*angle_fn2(lp),...
        init_rot(4)-cnst*angle_fn2(lp));        
    min_act(:,lp) = ff4(b_range(lp),init_rot(1)+minimum3(lp),...
        init_rot(2)-minimum3(lp),init_rot(3)+cnst*minimum3(lp),...
        init_rot(4)-cnst*minimum3(lp));         
end
%figure; plot(b_range,angle_fn)
%figure; plot(b_range,min_smooth)
figure; plot(b_range,1/30*min_smooth.'*[4,-1,-1;-1,4,-1;-1,-1,4])
hold on; plot(b_range,1/30*min_act.'*[4,-1,-1;-1,4,-1;-1,-1,4],'x')
%% New analysis
tmp = cat(1,fval2{:}); rnj = tmp < 2e-8;
shift2 = cat(1,minimum2{:});
angle_fn =zeros(length(b_range),2);  angle_fn2 = zeros(length(b_range),2);  
for lp=1:2
   % tmp2 = ezfit( b_range(rnj).',shift2(rnj,lp),'a*x^2+b*x^4'); 
   % angle_fn(:,lp) =  tmp2.m(1).*b_range.^2+tmp2.m(2).*b_range.^4;
   tmp2 = ezfit( b_range(rnj).',shift2(rnj,lp),'b*x^2'); 
    angle_fn(:,lp) =  tmp2.m(1).*b_range.^2;
    pp = splinefit(b_range(rnj),shift2(rnj,lp),3);
    angle_fn2(:,lp) = ppval(pp,b_range);
end
ff4 = matlabFunction([(e1.'*e2 .* e3.' *e4),(e1.'*e3 .* e2.' *e4),...
    e1.'*e4 .* e2.' *e3]);
min_smooth = zeros(3,length(b_range));
for lp = 1:length(b_range)

    min_smooth(:,lp) = ff4(b_range(lp),init_rot(1)+angle_fn2(lp,1),...
        init_rot(2)-angle_fn2(lp,1),init_rot(3)+angle_fn2(lp,2),...
        init_rot(4)-angle_fn2(lp,2));  
    min_act(:,lp) = ff4(b_range(lp),init_rot(1)+minimum2{lp}(1),...
        init_rot(2)-minimum2{lp}(1),init_rot(3)+minimum2{lp}(2),...
        init_rot(4)-minimum2{lp}(2));         
end
figure; plot(b_range,angle_fn2)
figure; plot1 = plot(b_range,min_smooth);

set(plot1(1),'DisplayName','$(\mathbf{p}_1 \cdot \mathbf{p}_2)  (\mathbf{p}_3 \cdot \mathbf{p}_4  )$');
set(plot1(2),'LineStyle','--',...
    'DisplayName','$(\mathbf{p}_1 \cdot \mathbf{p}_3)  (\mathbf{p}_2 \cdot \mathbf{p}_4  )$');
set(plot1(3),'LineStyle','-.',...
    'DisplayName','$(\mathbf{p}_1 \cdot \mathbf{p}_4)  (\mathbf{p}_3 \cdot \mathbf{p}_2  )$');

xlabel('\theta_2','FontSize',16);
ylabel('Amplitude','FontSize',16);
legend show

figure;plot2 =  plot(b_range,1/30*min_smooth.'*[4,-1,-1;-1,4,-1;-1,-1,4]);
xlabel('\theta_2','FontSize',16);
ylabel('Amplitude','FontSize',16);
set(plot2(1),'DisplayName','Comp 1');
set(plot2(2),'LineStyle','--','DisplayName','Comp 2');
set(plot2(3),'LineStyle','-.','DisplayName','Comp 3');
legend show
%hold on; plot(b_range,1/30*min_act.'*[4,-1,-1;-1,4,-1;-1,-1,4])

for lp = 1:length(b_range)
c = b_range(lp).^2/4; c2 = b_range(lp).^4/4/sqrt(2);
    testtt(lp,:)= f(b_range(lp),c+9*c2+init_rot(1),-c-9*c2+init_rot(2),...
                3*c+134*c2+init_rot(3),-3*c-134*c2+init_rot(4));
end

%% Analyse
shift2 = cat(1,minimum2{:});%-repmat(minimum2{1},length(b_range),1);
shift = cat(1,minimum{:});%-repmat(minimum{1},length(b_range),1);
order  = 2;
rnj = cat(1,fval2{:})< 3e-6 * max(0.5,b_range.'*100); %good points;
%figure;plot(b_range(rnj),[shift(rnj,:),shift2(rnj,:)])
shift_fit = zeros(4,order+1);
[angle_fn,angle_fn2,shift_smooth] =deal(zeros(length(b_range),4)); 
shift_fit2 = zeros(4,1);
b_long = pi*linspace(0,1/2,1000); shift_long =zeros(length(b_long),4);
for lp2= 1:4
   shift_fit(lp2,:) = polyfit( b_range(rnj).',shift2(rnj,lp2),order);
   tmp = ezfit( b_range(rnj).',shift2(rnj,lp2),'a*x^2');
   shift_fit2(lp2,:) =  tmp.m; 
   angle_fn(:,lp2) =  tmp.m.*b_range.^2;
   angle_fn2(:,lp2) = minimum{2}(lp2) + tmp.m.*b_range.^2;
   %shift_fit(lp2,1)=0; %cant just set this to zero..
   shift_smooth(:,lp2) = polyval(shift_fit(lp2,:),b_range.');
   shift_long(:,lp2) = polyval(shift_fit(lp2,:),b_long.');
end

min_smooth = b_range*0; min_smooth2 = zeros(3,length(b_range));%b_long*0;
fff =matlabFunction(abs(e1.'*e2 .* e3.' *e4) + abs(e1.'*e3 .* e2.' *e4)+...
    abs(e1.'*e4 .* e2.' *e3));
fff2 = matlabFunction([(e1.'*e2 .* e3.' *e4),(e1.'*e3 .* e2.' *e4),...
    e1.'*e4 .* e2.' *e3]);
for lp = 1:length(b_range)
    
    %min_smooth(lp) = fff(b_range(lp),minimum{1}(1)+shift_smooth(lp,1),minimum{1}(2)+shift_smooth(lp,2),...
    %minimum{1}(3)+shift_smooth(lp,3),minimum{1}(4)+shift_smooth(lp,4));
    min_smooth(lp) = fff(b_range(lp),angle_fn2(lp,1),angle_fn2(lp,2),...
                            angle_fn2(lp,3),angle_fn2(lp,4));
    min_smooth2(:,lp) = fff2(b_range(lp),angle_fn2(lp,1),angle_fn2(lp,2),...
                            angle_fn2(lp,3),angle_fn2(lp,4));                        
                        
end
% for lp = 1:length(b_long)
%     
%     min_smooth2(lp) = fff(b_long(lp),minimum{1}(1)+shift_long(lp,1),minimum{1}(2)+shift_long(lp,2),...
%     minimum{1}(3)+shift_long(lp,3),minimum{1}(4)+shift_long(lp,4));
% end

figure; plot(b_range,min_smooth2)
figure; plot(b_range,angle_fn)
%%

for lp = 1:length(b_range)
e1min = double(subs(e1,[t1,b],[minimum2{lp}(1),b_range(lp)]));
e2min = double(subs(e2,[t2,b],[minimum2{lp}(2),b_range(lp)]));
e3min = double(subs(e3,[t3,b],[minimum2{lp}(3),b_range(lp)]));
e4min = double(subs(e4,[t4,b],[minimum2{lp}(4),b_range(lp)]));
AA=[e1min,e2min,e3min,e4min]; 
AAx(lp,:) = AA(1,:); AAy(lp,:) = AA(2,:);AAz(lp,:) = AA(3,:);
end
%fit polynomials to these
%%
order=3;
[x_fit,y_fit,z_fit] = deal(zeros(4,order+1));
[AAx_fit,AAy_fit,AAz_fit] =  deal(zeros(length(b_range),4));
rnj = cat(1,fval2{:})<2e-6;%.*b_range.'*100 ; %keep only good solutions
for lp2= 1:4
   x_fit(lp2,:) = polyfit( b_range(rnj).',AAx(rnj,lp2),order);
   y_fit(lp2,:) = polyfit( b_range(rnj).',AAy(rnj,lp2),order); 
   z_fit(lp2,:) = polyfit( b_range(rnj).',AAz(rnj,lp2),order); 
   AAx_fit(:,lp2) = polyval(x_fit(lp2,:),b_range.');
   AAy_fit(:,lp2) = polyval(y_fit(lp2,:),b_range.');
   AAz_fit(:,lp2) = polyval(z_fit(lp2,:),b_range.');
end
AA_fit = cat(3,AAx_fit,AAy_fit,AAz_fit);
   fval3 = abs(sum(AA_fit(:,1,:).*AA_fit(:,2,:),3).* sum(AA_fit(:,3,:).*AA_fit(:,4,:),3))+...
            abs(sum(AA_fit(:,1,:).*AA_fit(:,3,:),3).* sum(AA_fit(:,2,:).*AA_fit(:,4,:),3))+...
            abs(sum(AA_fit(:,3,:).*AA_fit(:,2,:),3).* sum(AA_fit(:,1,:).*AA_fit(:,4,:),3));
figure; plot(b_range,[cat(1,fval{:}),cat(1,fval2{:}),fval3])
%%
 figure; subplot(3,1,1); plot(b_range(1:end),[AAx(1:end,1)-1,AAx(1:end,2:4)]); 
 ylabel('\Delta \epsilon_x')
 subplot(3,1,2); plot(b_range(1:end),[AAy(1:end,1),AAy(1:end,2:4)-1]);  
ylabel('\Delta \epsilon_y')
 subplot(3,1,3); plot(b_range(1:end),AAz(1:end,:)); 
  xlabel('Radial angle (radians)'); ylabel('\Delta \epsilon_z')
%plot(b_range(1:end/2),AAy(1:end/2,:),'--')
%plot(b_range(1:end/2),AAz(1:end/2,:),'.-')
%%
figure; subplot(3,1,1); plot(b_range(1:end),[AAx_fit(1:end,1)-1,AAx_fit(1:end,2:4)]); 
 ylabel('\Delta \epsilon_x')
subplot(3,1,2); plot(b_range(1:end),[AAy_fit(1:end,1),AAy_fit(1:end,2:4)-1]);  
ylabel('\Delta \epsilon_y')
subplot(3,1,3); plot(b_range(1:end),AAz_fit(1:end,:)); 
 xlabel('Radial angle (radians)'); ylabel('\Delta \epsilon_z')

 
 %% Alternative minimisation using polynomials
 
 xx = sym('x','real');
 
 syms a1 a2 a3 a4 b1 b2 b3 b4 c1 c2 c3 c4
 
 e1x = 1+a1*xx^2;  e2x = a2*xx^2;  e3x = a3*xx^2;  e4x = a4*xx^2;
 e1y = b1*xx^2;  e2y = 1 + b2*xx^2;  e3y = 1 + b3*xx^2;  e4y = 1 + b4*xx^2;
 e1z = c1*xx;  e2z = c2*xx;  e3z = c3*xx;  e4z = c4*xx;
 
 to_min = abs((e1x*e2x + e1y*e2y + e1z*e2z) * (e3x*e4x + e3y*e4y  + e3z*e4z )) +...
 abs((e1x*e3x + e1y*e3y + e1z*e3z) * (e2x*e4x + e2y*e4y  + e2z*e4z )) +...
 abs((e1x*e4x + e1y*e4y + e1z*e4z) * (e1x*e4x + e1y*e4y  + e1z*e4z )) ;
 
const = [1==e1x^2 + e1y^2 + e1z^2,1==e2x^2 + e2y^2 + e2z^2,...
        1==e3x^2 + e3y^2 + e3z^2, 1==e4x^2 + e4y^2 + e4z^2];

%%
    tset = [angle_fn(:,1),-angle_fn(:,1),angle_fn(:,2),-angle_fn(:,2)];
    tset = tset + repmat(init_rot,length(b_range),1);
    clear T T1 T2 T3 T4
for lp = 1:length(b_range)
e1min = double(subs(e1,[t1,b],[tset(lp,1),b_range(lp)]));
e2min = double(subs(e2,[t2,b],[tset(lp,2),b_range(lp)]));
e3min = double(subs(e3,[t3,b],[tset(lp,3),b_range(lp)]));
e4min = double(subs(e4,[t4,b],[tset(lp,4),b_range(lp)]));
AA=[e1min,e2min,e3min,e4min]; 
AAx(lp,:) = AA(1,:); AAy(lp,:) = AA(2,:);AAz(lp,:) = AA(3,:);
end
dAAx = AAx -repmat(AAx(1,:),length(b_range),1);
dAAy = AAy -repmat(AAy(1,:),length(b_range),1);
dAAz = AAz -repmat(AAz(1,:),length(b_range),1);
%%

tset = [angle_fn(:,1),-angle_fn(:,1),angle_fn(:,2),-angle_fn(:,2)];
clear T T1 T2 T3 T4
for lp= 1:size(tset,1)
Test1 = simplify(subs(e1.'*e2 .* e3.' *e4,[t1,t2,t3,t4],tset(lp,:)));
Test2 = simplify(subs(e1.'*e3 .* e2.' *e4,[t1,t2,t3,t4],tset(lp,:)));
Test3 = simplify(subs(e1.'*e4 .* e2.' *e3,[t1,t2,t3,t4],tset(lp,:)));
T(lp,:)=[Test1,Test2,Test3];

Test11 = simplify(subs(cross(k1,e1).'*e2 .* e3.' *e4,[t1,t2,t3,t4],tset(lp,:)));
Test21 = simplify(subs(cross(k1,e1).'*e3 .* e2.' *e4,[t1,t2,t3,t4],tset(lp,:)));
Test31 = simplify(subs(cross(k1,e1).'*e4 .* e2.' *e3,[t1,t2,t3,t4],tset(lp,:)));
T1(lp,:)=[Test11,Test21,Test31];

Test12 = simplify(subs(e1.'*cross(k2,e2) .* e3.' *e4,[t1,t2,t3,t4],tset(lp,:)));
Test22 = simplify(subs(e1.'*e3 .* cross(k2,e2).' *e4,[t1,t2,t3,t4],tset(lp,:)));
Test32 = simplify(subs(e1.'*e4 .* cross(k2,e2).' *e3,[t1,t2,t3,t4],tset(lp,:)));
T2(lp,:)=[Test12,Test22,Test32];

Test13 = simplify(subs(e1.'*e2 .* cross(k3,e3).' *e4,[t1,t2,t3,t4],tset(lp,:)));
Test23 = simplify(subs(e1.'*cross(k3,e3) .* e2.' *e4,[t1,t2,t3,t4],tset(lp,:)));
Test33 = simplify(subs(e1.'*e4 .* e2.'*cross(k3,e3),[t1,t2,t3,t4],tset(lp,:)));
T3(lp,:)=[Test13,Test23,Test33];

Test12 = simplify(subs(e1.'*e2 .* e3.' *cross(k4,e4),[t1,t2,t3,t4],tset(lp,:)));
Test22 = simplify(subs(e1.'*e3 .* e2.' *cross(k4,e4),[t1,t2,t3,t4],tset(lp,:)));
Test32 = simplify(subs(e1.'*cross(k4,e4) .* e2.' * e3,[t1,t2,t3,t4],tset(lp,:)));
T4(lp,:)=[Test13,Test23,Test33];
end


%% Plot graphs of all the angles between the polarizat
