syms a b c t1 t2 t3 t4
ry = [cos(b),0,sin(b);0,1,0;-sin(b),0,cos(b)];
rz = [cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
rz2 = [cos(c),-sin(c),0;sin(c),cos(c),0;0,0,1];

rot_op = rz2*ry*rz;

e1 = subs(rot_op*[1;0;0],[c,a],[pi/4,t1]);
e2 = subs(rot_op*[1;0;0],[c,a],[3*pi/4,t2]);
e3 = subs(rot_op*[1;0;0],[c,a],[7*pi/4,t3]);
e4 = subs(rot_op*[1;0;0],[c,a],[5*pi/4,t4]);

k1 = subs(rot_op*[0;0;1],[c,a],[pi/4,t1]);
k2 = subs(rot_op*[0;0;1],[c,a],[3*pi/4,t2]);
k3 = subs(rot_op*[0;0;1],[c,a],[7*pi/4,t3]);
k4 = subs(rot_op*[0;0;1],[c,a],[5*pi/4,t4]);

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
%%

b_range = pi*linspace(0,1/32,300);

OPTIONS.InitTemp = 1e-8; OPTIONS.CoolSched = @(T) 0.9.*T;
OPTIONS.StopTemp = 1e-10; OPTIONS.StopVal = 1e-8;
OPTIONS.MaxConsRej = 4000;  OPTIONS.MaxTries = 1000;
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
%%
ff4= matlabFunction([(e1.'*e2 .* e3.' *e4),(e1.'*e3 .* e2.' *e4),...
    e1.'*e4 .* e2.' *e3]);
%tmp = ezfit( b_range,minimum3,'a*x^2');
%angle_fn =  tmp.m.*b_range.^2;
tmp = ezfit( b_range.',minimum3,'a*x^2+b*x');
angle_fn =  tmp.m(1).*b_range.^2+tmp.m(2).*b_range;
min_smooth = zeros(3,length(b_range));
for lp = 1:length(b_range)
    min_smooth(:,lp) = ff3(b_range(lp),init_rot(1)+angle_fn(lp),...
        init_rot(2)+cnst*angle_fn(lp),init_rot(3)-angle_fn(lp),...
        init_rot(4)-cnst*angle_fn(lp));                                        
end
figure; plot(b_range,angle_fn)
figure; plot(b_range,min_smooth)
%% New analysis
tmp = cat(1,fval2{:}); rnj = tmp < 1e-6;
shift2 = cat(1,minimum2{:});
angle_fn =zeros(length(b_range),2); 
for lp=1:2
    tmp2 = ezfit( b_range(rnj).',shift2(rnj,lp),'a*x^2'); 
    angle_fn(:,lp) =  tmp2.m.*b_range.^2;
end
ff3 = matlabFunction([(e1.'*e2 .* e3.' *e4),(e1.'*e3 .* e2.' *e4),...
    e1.'*e4 .* e2.' *e3]);
min_smooth = zeros(3,length(b_range));
for lp = 1:length(b_range)

    min_smooth(:,lp) = ff3(b_range(lp),init_rot(1)+angle_fn(lp,1),...
        init_rot(2)+angle_fn(lp,2),init_rot(3)-angle_fn(lp,1),...
        init_rot(4)-angle_fn(lp,2));                                        
end
figure; plot(b_range,min_smooth)

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

tset = minimum{lp};
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