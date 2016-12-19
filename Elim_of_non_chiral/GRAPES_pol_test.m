syms a b c t1 t2 t3 t4 b1 b2 cc real
ry = [cos(b),0,sin(b);0,1,0;-sin(b),0,cos(b)];
rz = [cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
rz2 = [cos(c),-sin(c),0;sin(c),cos(c),0;0,0,1];

rot_op = rz2*ry*rz;

%Contruct Grapes k vectors, with free parameters for e1,..,e4

e0 = rot_op*[1;0;0]; k0 = rot_op*[0;0;1];

e1 = subs(e0,[c,b,a],[-cc,b1,t1]); %cc == phi = acos(b2/b1)
e2 = subs(e0,[c,b,a],[0,b2,t2]);
e3 = subs(e0,[c,b,a],[pi,b2,t3]);
e4 = subs(e0,[c,b,a],[pi-cc,b1,t4]);

rot_inv = subs(rot_op,[c,b,a],[-t4,-b1,-pi+cc]); %rotates 4th back to x

k1 = subs(k0,[c,b,a],[-cc,b1,t1]);
k2 = subs(k0,[c,b,a],[0,b2,t2]); %k2 = subs(k0,[c,b,a],[0,b1*cos(cc),t2]);
k3 = subs(k0,[c,b,a],[pi,b2,t3]);%k3 = subs(k0,[c,b,a],[pi,b1*cos(cc),t3]);
k4 = subs(k0,[c,b,a],[pi-cc,b1,t4]);

%% t4 = 0 solutions
ee1 = subs(e1,t1,pi/2+cc + t1*b1^2); 
 ee2 = subs(e2,[b2,t2],[b1*cos(cc),t2*b1^2]); 
ee3 = subs(e3,[b2,t3],[b1*cos(cc),t3*b1^2]); 
 ee4 = subs(e4,t4,pi+cc + 0); %set t4 = 0
 
test = [ee1.'*ee2 * ee4.'*ee3,ee1.'*ee3 * ee2.'*ee4, ee1.'*ee4 * ee2.'*ee3];
test= test*[4,-1,-1;-1,4,-1;-1,-1,4]; test = test.'/30;

test2 = simplify(taylor(test,b1,0,'order',3));
test3 = solve(test2==[0;0;0],t1,t2,t3);

ee1 = subs(ee1,t1,test3.t1 + t1*b1^4); ee2 = subs(ee2,t2,test3.t2+ t2*b1^4); 
ee3 = subs(ee3,t3,test3.t3+ t3*b1^4); 
 
test = [ee1.'*ee2 * ee4.'*ee3,ee1.'*ee3 * ee2.'*ee4, ee1.'*ee4 * ee2.'*ee3];
test= test*[4,-1,-1;-1,4,-1;-1,-1,4]; test = test.'/30;

test2 = simplify(taylor(test,b1,0,'order',7));
test3 = solve(test2==[0;0;0],t1,t2,t3);

ee1 = subs(ee1,t1,test3.t1); ee2 = subs(ee2,t2,test3.t2); 
ee3 = subs(ee3,t3,test3.t3); 

test = [ee1.'*ee2 * ee4.'*ee3,ee1.'*ee3 * ee2.'*ee4, ee1.'*ee4 * ee2.'*ee3];
test= test*[4,-1,-1;-1,4,-1;-1,-1,4]; test = test.'/30;

simplify(taylor(test,b1,0,'order',7))
%%

%  init_rot = [pi/2+phi,0,pi,pi+phi]; 
%A = sym('A'); A2 = sym('B');
soc = sin(2*cc)*b1^2/8; syms D
%% Minimise deformation solutions
% ee1 = subs(e1,t1,pi/2+cc + 3/8*sin(2*cc)*b1^2); 
% ee2 = subs(e2,[b2,t2],[b1*cos(cc),5/8*sin(2*cc)*b1^2]); 
% ee3 = subs(e3,[b2,t3],[b1*cos(cc),pi-3/8*sin(2*cc)*b1^2]); 
% ee4 = subs(e4,t4,pi+cc - 5/8*sin(2*cc)*b1^2);
ee1 = subs(e1,t1,pi/2+cc + (D + 8)*soc ); 
ee2 = subs(e2,[b2,t2],[b1*cos(cc),+ (D + 10)*soc]); 
ee3 = subs(e3,[b2,t3],[b1*cos(cc),pi+ (D + 2)*soc]); 
ee4 = subs(e4,t4,pi+cc + D*soc);
test = [ee1.'*ee2 * ee4.'*ee3,ee1.'*ee3 * ee2.'*ee4 ,ee1.'*ee4 * ee2.'*ee3];
test= test*[4,-1,-1;-1,4,-1;-1,-1,4]; test = test.';

simplify(taylor(test,b1,0,'order',3))
simplify(taylor(test,b1,0,'order',5))

ee1 = subs(e1,t1,pi/2+cc + (-5 + 8)*soc +t1*b1^4); 
ee2 = subs(e2,[b2,t2],[b1*cos(cc),+ (-5 + 10)*soc+t2*b1^4]); 
ee3 = subs(e3,[b2,t3],[b1*cos(cc),pi+ (-5 + 2)*soc+t3*b1^4]); 
ee4 = subs(e4,t4,pi+cc + -5*soc+t4*b1^4);
test = [ee1.'*ee2 * ee4.'*ee3,ee1.'*ee3 * ee2.'*ee4 ,ee1.'*ee4 * ee2.'*ee3];
test= test*[4,-1,-1;-1,4,-1;-1,-1,4]; test = test.'/30;

test2 =simplify(taylor(test,b1,0,'order',5));
test3 = solve(test2 ==[0;0;0],t1,t2,t3);
test4=[test3.t1;test3.t2;test3.t3;t4];  %t4 = (37*sin(4*cc))/384 - (25*sin(2*cc))/64
test5 = subs(test4,t4,solve(sum(test4),t4));

ee12 = subs(ee1,t1,test5(1)); ee22 = subs(ee2,t2,test5(2));
ee32 = subs(ee3,t3,test5(3)); ee42 = subs(ee4,t4,test5(4));
test = [ee12.'*ee22 * ee42.'*ee32,ee12.'*ee32 * ee22.'*ee42, ...
        ee12.'*ee42 * ee22.'*ee32];
test= test*[4,-1,-1;-1,4,-1;-1,-1,4]; test = test.'/30;
simplify(taylor(test,b1,0,'order',7))
%%

Ee1 = subs(e1,t1,cc + D*soc); 
Ee2 = subs(e2,[b2,t2],[b1*cos(cc),pi/2 + (D+2)*soc]); 
Ee3 = subs(e3,[b2,t3],[b1*cos(cc),pi + (D+2)*soc]); 
Ee4 = subs(e4,t4,pi + cc + D*soc);
test2 = [Ee1.'*Ee2 * Ee4.'*Ee3,Ee1.'*Ee3 * Ee2.'*Ee4 Ee1.'*Ee4 * Ee2.'*Ee3];
test2= test2*[4,-1,-1;-1,4,-1;-1,-1,4]; test2 = test2.';

simplify(taylor(test2,b1,0,'order',3))
simplify(taylor(test2,b1,0,'order',5))
%% Order 2 min
EE1 = subs(e1,t1,cc+ D*soc); 
EE2 = subs(e2,[b2,t2],[b1*cos(cc),pi+ (D+2)*soc]); 
EE3 = subs(e3,[b2,t3],[b1*cos(cc),3*pi/2+ (D+2)*soc]);
EE4 = subs(e4,t4,pi+cc+ D*soc);
test3 = [EE1.'*EE2 * EE4.'*EE3,EE1.'*EE3 * EE2.'*EE4 EE1.'*EE4 * EE2.'*EE3];
test3= test3*[4,-1,-1;-1,4,-1;-1,-1,4]; test3 = test3.';

simplify(taylor(test3,b1,0,'order',3))
simplify(taylor(test3,b1,0,'order',5))

%% Order 4 min

ee1 = subs(e1,t1,pi/2+cc + (-5 + 8)*soc +t1*b1^4 ); 
ee2 = subs(e2,b2,b1*cos(cc)); ee2 =subs(ee2, t2,(-5 + 10)*soc +t2*b1^4 );
ee3 = subs(e3,b2,b1*cos(cc)); ee3 =subs(ee3,t3,pi+ (-5 + 2)*soc +t3*b1^4); 
ee4 = subs(e4,t4,pi+cc + -5*soc +t4*b1^4);
tmp= [ee1.'*ee2 * ee4.'*ee3,ee1.'*ee3 * ee2.'*ee4 ee1.'*ee4 * ee2.'*ee3];
tmp = tmp*[4,-1,-1;-1,4,-1;-1,-1,4]/30;
tmp2 = simplify(taylor(tmp,b1,0,'order',5));
tmp3 = solve(tmp2 ==[0,0,0],t1,t2,t3);
t4_sol = solve(t4 + tmp3.t1 + tmp3.t2 + tmp3.t3==0,t4);

ee1 = subs(ee1,t1,tmp3.t1); ee2 = subs(ee2,t2,tmp3.t2); 
ee3 = subs(ee3,t3,tmp3.t3);   
%test= [ee1.'*ee2 * ee4.'*ee3,ee1.'*ee3 * ee2.'*ee4 ee1.'*ee4 * ee2.'*ee3];
%test = test*[4,-1,-1;-1,4,-1;-1,-1,4]/30;
%simplify(taylor(test,b1,0,'order',3))
%simplify(taylor(test,b1,0,'order',5))
theta_fn = acos(subs([dot(ee1,ee2);dot(ee1,ee3);dot(ee1,ee4);...
                dot(ee2,ee3);dot(ee2,ee4);dot(ee3,ee4)],t4,t4_sol));
Theta_fn = matlabFunction(theta_fn);

%%
cnt = 0; %test impact of random noise
namp_rng = 1/60; namp_rng = repmat(namp_rng,[1,1000]);
for lp = 1:length(namp_rng)

namp = namp_rng(lp)*pi/180; R = namp*(rand(4)-1/2); %random noise

ee1n = subs(e1,t1,pi/2+cc + (-5 + 8)*soc +tmp3.t1*b1^4 +R(1)); 
ee2n = subs(e2,b2,b1*cos(cc)); ee2n =subs(ee2n, t2,(-5 + 10)*soc +tmp3.t2*b1^4+R(2));
ee3n = subs(e3,b2,b1*cos(cc)); ee3n =subs(ee3n,t3,pi+ (-5 + 2)*soc +tmp3.t3*b1^4+R(3)); 
ee4n = subs(e4,t4,pi+cc + -5*soc +t4_sol*b1^4 + R(4)); 

%{
T1=-(16*sin(b)^4 - 8*sin(b)^2 - 22*b1^2*sin(b)^2 + 12*b1^2*sin(b)^4 - 8*b1^4*t4 + 8*b1^2 - 4*b1^2*sin(b)^2*sin(cc)^2 + 8*b1^2*sin(b)^4*sin(cc)^2 + 17*b1^4*sin(b)^2*sin(cc)^2 - 17*b1^4*sin(b)^2*sin(cc)^4 - 19*b1^4*sin(b)^4*sin(cc)^2 + 19*b1^4*sin(b)^4*sin(cc)^4 + 16*b1^4*t4*sin(b)^2 + 16*b1^4*t4*sin(b)^2*sin(cc)^2 - 32*b1^4*t4*sin(b)^4*sin(cc)^2)/(8*b1^4*cos(2*b)*(2*sin(b)^2*sin(cc)^2 - 2*sin(b)^2 + 1));
T2=                                  -(16*cos(b)^2 - 48*cos(b)^4 + 32*cos(b)^6 - 32*cos(cc)^2 + 32*cos(cc)^4 - 40*b1^2*cos(cc)^2 + 120*b1^2*cos(cc)^4 - 80*b1^2*cos(cc)^6 + 38*b1^4*cos(cc)^4 - 136*b1^4*cos(cc)^6 + 158*b1^4*cos(cc)^8 - 60*b1^4*cos(cc)^10 + 96*cos(b)^2*cos(cc)^2 - 128*cos(b)^2*cos(cc)^4 - 64*cos(b)^4*cos(cc)^2 + 192*cos(b)^4*cos(cc)^4 - 160*cos(b)^6*cos(cc)^4 + 64*cos(b)^8*cos(cc)^4 + 48*cos(b)^3*cos(cc)*sin(cc) - 16*cos(b)^5*cos(cc)*sin(cc) - 32*cos(b)^7*cos(cc)*sin(cc) + 160*b1^2*cos(b)^2*cos(cc)^2 - 480*b1^2*cos(b)^2*cos(cc)^4 - 240*b1^2*cos(b)^4*cos(cc)^2 - 29*b1^4*cos(b)^2*cos(cc)^2 + 320*b1^2*cos(b)^2*cos(cc)^6 + 720*b1^2*cos(b)^4*cos(cc)^4 + 200*b1^2*cos(b)^6*cos(cc)^2 - 19*b1^4*cos(b)^2*cos(cc)^4 + 87*b1^4*cos(b)^4*cos(cc)^2 - 480*b1^2*cos(b)^4*cos(cc)^6 - 600*b1^2*cos(b)^6*cos(cc)^4 - 80*b1^2*cos(b)^8*cos(cc)^2 + 258*b1^4*cos(b)^2*cos(cc)^6 - 229*b1^4*cos(b)^4*cos(cc)^4 - 58*b1^4*cos(b)^6*cos(cc)^2 + 400*b1^2*cos(b)^6*cos(cc)^6 + 240*b1^2*cos(b)^8*cos(cc)^4 - 330*b1^4*cos(b)^2*cos(cc)^8 + 132*b1^4*cos(b)^4*cos(cc)^6 + 250*b1^4*cos(b)^6*cos(cc)^4 - 160*b1^2*cos(b)^8*cos(cc)^6 + 120*b1^4*cos(b)^2*cos(cc)^10 + 10*b1^4*cos(b)^4*cos(cc)^8 - 342*b1^4*cos(b)^6*cos(cc)^6 - 40*b1^4*cos(b)^8*cos(cc)^4 + 270*b1^4*cos(b)^6*cos(cc)^8 + 88*b1^4*cos(b)^8*cos(cc)^6 - 120*b1^4*cos(b)^6*cos(cc)^10 - 108*b1^4*cos(b)^8*cos(cc)^8 + 60*b1^4*cos(b)^8*cos(cc)^10 + 64*b1^4*t4*cos(cc)^2 - 192*b1^4*t4*cos(cc)^4 + 128*b1^4*t4*cos(cc)^6 - 16*cos(b)*cos(cc)*sin(cc) - 200*b1^2*cos(b)^3*cos(cc)^3*sin(cc) + 80*b1^2*cos(b)^3*cos(cc)^5*sin(cc) + 240*b1^2*cos(b)^5*cos(cc)^3*sin(cc) - 76*b1^4*cos(b)^3*cos(cc)^3*sin(cc) - 240*b1^2*cos(b)^5*cos(cc)^5*sin(cc) + 174*b1^4*cos(b)^3*cos(cc)^5*sin(cc) + 27*b1^4*cos(b)^5*cos(cc)^3*sin(cc) + 160*b1^2*cos(b)^7*cos(cc)^5*sin(cc) - 98*b1^4*cos(b)^3*cos(cc)^7*sin(cc) - 81*b1^4*cos(b)^5*cos(cc)^5*sin(cc) + 38*b1^4*cos(b)^7*cos(cc)^3*sin(cc) + 54*b1^4*cos(b)^5*cos(cc)^7*sin(cc) - 52*b1^4*cos(b)^7*cos(cc)^5*sin(cc) + 14*b1^4*cos(b)^7*cos(cc)^7*sin(cc) - 20*b1^2*cos(b)*cos(cc)*sin(cc) - 256*b1^4*t4*cos(b)^2*cos(cc)^2 + 768*b1^4*t4*cos(b)^2*cos(cc)^4 + 384*b1^4*t4*cos(b)^4*cos(cc)^2 - 512*b1^4*t4*cos(b)^2*cos(cc)^6 - 1152*b1^4*t4*cos(b)^4*cos(cc)^4 - 320*b1^4*t4*cos(b)^6*cos(cc)^2 + 768*b1^4*t4*cos(b)^4*cos(cc)^6 + 960*b1^4*t4*cos(b)^6*cos(cc)^4 + 128*b1^4*t4*cos(b)^8*cos(cc)^2 - 640*b1^4*t4*cos(b)^6*cos(cc)^6 - 384*b1^4*t4*cos(b)^8*cos(cc)^4 + 256*b1^4*t4*cos(b)^8*cos(cc)^6 + 40*b1^2*cos(b)*cos(cc)^3*sin(cc) + 80*b1^2*cos(b)^3*cos(cc)*sin(cc) - 60*b1^2*cos(b)^5*cos(cc)*sin(cc) + 19*b1^4*cos(b)*cos(cc)^3*sin(cc) - 40*b1^2*cos(b)^7*cos(cc)*sin(cc) - 49*b1^4*cos(b)*cos(cc)^5*sin(cc) + 30*b1^4*cos(b)*cos(cc)^7*sin(cc) - 64*b1^4*t4*cos(b)*cos(cc)^3*sin(cc) - 96*b1^4*t4*cos(b)^3*cos(cc)*sin(cc) + 32*b1^4*t4*cos(b)^5*cos(cc)*sin(cc) + 64*b1^4*t4*cos(b)^7*cos(cc)*sin(cc) + 192*b1^4*t4*cos(b)^3*cos(cc)^3*sin(cc) - 64*b1^4*t4*cos(b)^5*cos(cc)^3*sin(cc) - 128*b1^4*t4*cos(b)^7*cos(cc)^3*sin(cc) + 32*b1^4*t4*cos(b)*cos(cc)*sin(cc))/(32*b1^4*cos(b)^3*cos(cc)*sin(cc)*(2*cos(b)^2 - 1)*(2*cos(cc)^2 - 1)*(2*cos(b)^2*cos(cc)^2 - 2*cos(cc)^2 + 1));
T3= (16*cos(b)^2 - 48*cos(b)^4 + 32*cos(b)^6 - 32*cos(cc)^2 + 32*cos(cc)^4 - 40*b1^2*cos(cc)^2 + 120*b1^2*cos(cc)^4 - 80*b1^2*cos(cc)^6 + 38*b1^4*cos(cc)^4 - 136*b1^4*cos(cc)^6 + 158*b1^4*cos(cc)^8 - 60*b1^4*cos(cc)^10 + 96*cos(b)^2*cos(cc)^2 - 128*cos(b)^2*cos(cc)^4 - 64*cos(b)^4*cos(cc)^2 + 192*cos(b)^4*cos(cc)^4 - 160*cos(b)^6*cos(cc)^4 + 64*cos(b)^8*cos(cc)^4 - 48*cos(b)^3*cos(cc)*sin(cc) + 16*cos(b)^5*cos(cc)*sin(cc) + 32*cos(b)^7*cos(cc)*sin(cc) + 160*b1^2*cos(b)^2*cos(cc)^2 - 480*b1^2*cos(b)^2*cos(cc)^4 - 240*b1^2*cos(b)^4*cos(cc)^2 - 29*b1^4*cos(b)^2*cos(cc)^2 + 320*b1^2*cos(b)^2*cos(cc)^6 + 720*b1^2*cos(b)^4*cos(cc)^4 + 200*b1^2*cos(b)^6*cos(cc)^2 - 19*b1^4*cos(b)^2*cos(cc)^4 + 87*b1^4*cos(b)^4*cos(cc)^2 - 480*b1^2*cos(b)^4*cos(cc)^6 - 600*b1^2*cos(b)^6*cos(cc)^4 - 80*b1^2*cos(b)^8*cos(cc)^2 + 258*b1^4*cos(b)^2*cos(cc)^6 - 229*b1^4*cos(b)^4*cos(cc)^4 - 58*b1^4*cos(b)^6*cos(cc)^2 + 400*b1^2*cos(b)^6*cos(cc)^6 + 240*b1^2*cos(b)^8*cos(cc)^4 - 330*b1^4*cos(b)^2*cos(cc)^8 + 132*b1^4*cos(b)^4*cos(cc)^6 + 250*b1^4*cos(b)^6*cos(cc)^4 - 160*b1^2*cos(b)^8*cos(cc)^6 + 120*b1^4*cos(b)^2*cos(cc)^10 + 10*b1^4*cos(b)^4*cos(cc)^8 - 342*b1^4*cos(b)^6*cos(cc)^6 - 40*b1^4*cos(b)^8*cos(cc)^4 + 270*b1^4*cos(b)^6*cos(cc)^8 + 88*b1^4*cos(b)^8*cos(cc)^6 - 120*b1^4*cos(b)^6*cos(cc)^10 - 108*b1^4*cos(b)^8*cos(cc)^8 + 60*b1^4*cos(b)^8*cos(cc)^10 + 64*b1^4*t4*cos(cc)^2 - 192*b1^4*t4*cos(cc)^4 + 128*b1^4*t4*cos(cc)^6 + 16*cos(b)*cos(cc)*sin(cc) + 72*b1^2*cos(b)^3*cos(cc)^3*sin(cc) + 48*b1^2*cos(b)^3*cos(cc)^5*sin(cc) + 80*b1^2*cos(b)^5*cos(cc)^3*sin(cc) + 76*b1^4*cos(b)^3*cos(cc)^3*sin(cc) - 144*b1^2*cos(b)^5*cos(cc)^5*sin(cc) - 128*b1^2*cos(b)^7*cos(cc)^3*sin(cc) - 174*b1^4*cos(b)^3*cos(cc)^5*sin(cc) - 27*b1^4*cos(b)^5*cos(cc)^3*sin(cc) + 96*b1^2*cos(b)^7*cos(cc)^5*sin(cc) + 98*b1^4*cos(b)^3*cos(cc)^7*sin(cc) + 81*b1^4*cos(b)^5*cos(cc)^5*sin(cc) - 38*b1^4*cos(b)^7*cos(cc)^3*sin(cc) - 54*b1^4*cos(b)^5*cos(cc)^7*sin(cc) + 52*b1^4*cos(b)^7*cos(cc)^5*sin(cc) - 14*b1^4*cos(b)^7*cos(cc)^7*sin(cc) + 20*b1^2*cos(b)*cos(cc)*sin(cc) - 256*b1^4*t4*cos(b)^2*cos(cc)^2 + 768*b1^4*t4*cos(b)^2*cos(cc)^4 + 384*b1^4*t4*cos(b)^4*cos(cc)^2 - 512*b1^4*t4*cos(b)^2*cos(cc)^6 - 1152*b1^4*t4*cos(b)^4*cos(cc)^4 - 320*b1^4*t4*cos(b)^6*cos(cc)^2 + 768*b1^4*t4*cos(b)^4*cos(cc)^6 + 960*b1^4*t4*cos(b)^6*cos(cc)^4 + 128*b1^4*t4*cos(b)^8*cos(cc)^2 - 640*b1^4*t4*cos(b)^6*cos(cc)^6 - 384*b1^4*t4*cos(b)^8*cos(cc)^4 + 256*b1^4*t4*cos(b)^8*cos(cc)^6 - 40*b1^2*cos(b)*cos(cc)^3*sin(cc) - 48*b1^2*cos(b)^3*cos(cc)*sin(cc) - 4*b1^2*cos(b)^5*cos(cc)*sin(cc) - 19*b1^4*cos(b)*cos(cc)^3*sin(cc) + 40*b1^2*cos(b)^7*cos(cc)*sin(cc) + 49*b1^4*cos(b)*cos(cc)^5*sin(cc) - 30*b1^4*cos(b)*cos(cc)^7*sin(cc) + 64*b1^4*t4*cos(b)*cos(cc)^3*sin(cc) + 96*b1^4*t4*cos(b)^3*cos(cc)*sin(cc) - 32*b1^4*t4*cos(b)^5*cos(cc)*sin(cc) - 64*b1^4*t4*cos(b)^7*cos(cc)*sin(cc) - 192*b1^4*t4*cos(b)^3*cos(cc)^3*sin(cc) + 64*b1^4*t4*cos(b)^5*cos(cc)^3*sin(cc) + 128*b1^4*t4*cos(b)^7*cos(cc)^3*sin(cc) - 32*b1^4*t4*cos(b)*cos(cc)*sin(cc))/(32*b1^4*cos(b)^3*cos(cc)*sin(cc)*(2*cos(b)^2 - 1)*(2*cos(cc)^2 - 1)*(2*cos(b)^2*cos(cc)^2 - 2*cos(cc)^2 + 1));
tmp2 = subs(tmp,[t1,t2,t3],[T1,T2,T3]);
simplify(taylor(tmp2,b1,0,'order',5))
%}

b_range = pi*linspace(0,1/20,3);

%tmp= [ee1.'*ee2 * ee4.'*ee3,ee1.'*ee3 * ee2.'*ee4 ee1.'*ee4 * ee2.'*ee3];
tmp2= [ee1n.'*ee2n * ee4n.'*ee3n,ee1n.'*ee3n * ee2n.'*ee4n ee1n.'*ee4n * ee2n.'*ee3n];
%tmp = subs(tmp,t4,t4_sol)*[4,-1,-1;-1,4,-1;-1,-1,4]/30;
tmp_noise = subs(tmp2,t4,t4_sol)*[4,-1,-1;-1,4,-1;-1,-1,4]/30;

%ff = matlabFunction(sum(abs(tmp)));
ff2 = matlabFunction(sum(abs(tmp_noise)));

inter_angle = zeros(length(b_range),6);
 for k = 1:length(b_range);
     b = b_range(k); phi = 0.23*pi;
   %  inter_angle(k,:) = Theta_fn(b,phi);
   %  error_in_cancel(k) = ff(b,phi);
     error_w_noise(k) = ff2(b,phi);
 end   
 cnt = cnt+1;
 error_w_noise_save(lp,:) = error_w_noise;
end
 %% Plot angles between Polarizations of the beams
 figure; plotyy(b_range/pi*180,(inter_angle(:,[1:3,6])-repmat(...
     inter_angle(1,[1:3,6]),length(b_range),1))/pi*180,...
     b_range/pi*180,(inter_angle(:,4:5)-repmat(inter_angle(1,4:5),length(b_range),1))/pi*180)
 
 %%
 figure1 = figure;  axes1 = axes('Parent',figure1);  hold(axes1,'on');
% plot1 = plot(b_range/pi*180,...
%     asinh(1000*(inter_angle-repmat(inter_angle(1,:),length(b_range),1))),'Parent',axes1);
% tmp = get(gca, 'YTick'); set (gca, 'YTickLabel', num2str(sinh(tmp(:))/1000/pi*180 , '%g')) 
% h = zoom; %define handle for 'zoom'  
% %action to be called right after zooming  
% set(h,'ActionPostCallback', {@mypostcallback}); 
plot1 = plot(b_range/pi*180,(inter_angle-repmat(inter_angle(1,:),length(b_range),1))/pi*180,'Parent',axes1);
set(plot1(1),'DisplayName','\theta_{12}');
set(plot1(2),'DisplayName','\theta_{13}');
set(plot1(3),'DisplayName','\theta_{14}','LineStyle','--');
set(plot1(4),'DisplayName','\theta_{23}');
set(plot1(5),'DisplayName','\theta_{24}','LineStyle',':');
set(plot1(6),'DisplayName','\theta_{34}','LineStyle','-.');

xlabel('\Theta (degrees)');  % xlim(axes1,[0 0.5]);
ylabel('Angle between polarization (degrees)');

box(axes1,'on');
legend1 = legend(axes1,'show');    set(legend1,'Position',...
[0.138834329030446 0.485804414711539 0.17092337735978 0.424815971329887]); 

figure2 = figure; axes2= axes('Parent',figure2);  hold(axes2,'on'); 
plot2 = plot(b_range/pi*180,(inter_angle(:,1:6)-repmat(inter_angle(1,1:6),length(b_range),1))/pi*180,'Parent',axes2);
[h_main, h_inset]=inset_figure(figure1, figure2,0.4);

%figure; plotyy(b_range,inter_angle-repmat(inter_angle(1,:),length(b_range),1),...
%        b_range,error_in_cancel);

%%
r_range = pi*linspace(0,1/32,50);
phi_range = -pi*linspace(0.01,1.6,10)/4;
minimum = zeros(4,length(phi_range),length(r_range)); 
fval = zeros(length(phi_range),length(r_range));

options.TolFun = 1e-8;  options.TolX = 1e-8; options.verbose=0;

ee1 = subs(e1,t1,pi/2+cc + (-5 + 8)*soc); 
ee2 = subs(e2,b2,b1*cos(cc)); ee2 =subs(ee2, t2,(-5 + 10)*soc);
ee3 = subs(e3,b2,b1*cos(cc)); ee3 =subs(ee3,t3,pi+ (-5 + 2)*soc ); 
ee4 = subs(e4,t4,pi+cc + -5*soc);
tmp= [ee1.'*ee2 * ee4.'*ee3,ee1.'*ee3 * ee2.'*ee4 ee1.'*ee4 * ee2.'*ee3];
tmp = tmp*[4,-1,-1;-1,4,-1;-1,-1,4]/30;

ff = matlabFunction(sum(abs(tmp)));
%{
for lp1 = 1:length(phi_range)
    
            phi = phi_range(lp1); %angle between beams 1 and 2
 for lp2 = 1:length(r_range)
     rr = r_range(lp2);
     ff2 = @(tt) ff(rr,phi,tt(1),tt(2),tt(3),tt(4));
    if lp2 > 1; guess = minimum(:,lp1,lp2-1); 
    else; guess = [0,0,0,0];   end
     [minimum(:,lp1,lp2),fval(lp1,lp2)] =  fminunc(ff2 ,guess,options);  
 end
 
end
%%
tmp = ezfit(squeeze(minimum(1,5,:)),'a*x^4+b*x^6');
%}
%%
%wtf was this section meant to do!?
r_range = pi*linspace(0,1/32,50);
phi_range = pi*linspace(0.01,1.6,10)/4;

ee_save = zeros(3,4,length(phi_range));
ff_save = zeros(3,length(phi_range));
for lp1 = 1:length(phi_range)
    rr = r_range(25);
            phi = phi_range(lp1); %angle between beams 1 and 2
 %init_rot = [(sin(phi) + 1)/cos(phi),0,pi,-(cos(phi) + 1)/sin(phi)]; 
 init_rot = [pi/2+phi,0,pi,pi+phi]; 
 e10 = subs(e1,[cc,b1,t1],[phi,rr,init_rot(1)]);
 e20 = subs(e2,[cc,b2,t2],[phi,rr*cos(phi),init_rot(2)]);
 e30 = subs(e3,[cc,b2,t3],[phi,rr*cos(phi),init_rot(3)]);
 e40 = subs(e4,[cc,b1,t4],[phi,rr,init_rot(4)]);
 ee_save(:,:,lp1) = double([e10,e20,e30,e40]);
 test2 =  subs(test,[cc],[phi]);
 for lp2 = 1:length(r_range)
test3 = subs(test2,[b1],[r_range(lp2)]);
ff_save(:,lp1,lp2) = test3;
 end
 ee_save(:,:,lp1) = double([e10,e20,e30,e40]);
end
ff_save2 = squeeze(sum(ff_save))/30;




%%
%{
OPTIONS.InitTemp = 1; OPTIONS.CoolSched = @(T) 0.9.*T;
OPTIONS.StopTemp = 1e-10; OPTIONS.StopVal = 1e-8;
OPTIONS.MaxConsRej = 4000;  OPTIONS.MaxTries = 1000;
for lp1 = 1:length(r_range)
    for lp2 = 1:length(phi_range)
        
        r_small = r_range(lp1)*cos(phi_range(lp2)); %smaller one
    
    tmp = subs(abs(e1.'*e2 .* e3.' *e4) + abs(e1.'*e3 .* e2.' *e4)+...
    abs(e1.'*e4 .* e2.' *e3),'b',b_range(lp));
f = matlabFunction(tmp);
ff = @(x) f(x(1),x(2),x(3),x(4));
    if lp > 1; guess = minimum{lp-1};
    else; guess = pi/4*[-1,-1,3,-3]; end %pi/4*[1,-1,3,3];
[minimum{lp},fval{lp}] = anneal(ff,guess,OPTIONS);
OPTIONS.InitTemp = OPTIONS.InitTemp*2;
    end
end
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

minimum = zeros(4,length(phi_range),length(r_range));
minimum2 = zeros(2,length(phi_range),length(r_range));
fval = zeros(length(phi_range),length(r_range)); fval2 = fval;

options.TolFun = 1e-8;  options.TolX = 1e-8; 

     tmp = [(e1.'*e2 .* e3.' *e4) , (e1.'*e3 .* e2.' *e4),...
            (e1.'*e4 .* e2.' *e3)];
    tmp2 = 1/30*tmp*[4,-1,-1;-1,4,-1;-1,-1,4] ; tmp3 = sum(abs(tmp2));
    f = matlabFunction(tmp3); %b1 b2 cc t1 t2 t3 t4
for lp1 = 1:length(phi_range)
    
            phi = phi_range(lp1); %angle between beams 1 and 2
            
  init_rot = [pi/2+phi,0,pi,pi+phi]; 
 
    for lp2 = 1:length(r_range)

  %init_rot = [acos(1/cos(r_range(lp2))),0,pi,-pi/2];   
        
        r_small = r_range(lp2)*cos(phi); %smaller one

ff = @(x) f(r_range(lp2),r_small,phi,...
    x(1)+init_rot(1),x(2)+init_rot(2),x(3)+init_rot(3),x(4)+init_rot(4));
ff2 = @(x) f(r_range(lp2),r_small,phi,...
    x(1)+init_rot(1),x(2)+init_rot(2),-x(1)+init_rot(3),-x(2)+init_rot(4));

    if lp2 > 1; guess = minimum(:,lp1,lp2-1).'; guess2 = minimum2(:,lp1,lp2-1).';
    else; guess = [0,0,0,0]; guess2 = [0,0];   end
%[minimum(:,lp1,lp2),fval(lp1,lp2)] =  fminunc(ff ,guess);    %more general
[minimum2(:,lp1,lp2),fval2(lp1,lp2)] =  fminunc(ff2 ,guess2);    %less general
    end
end

%% Analyse

min_smooth = fval*0; min_smooth2 = zeros([3,size(fval)]);%b_long*0;
fff = matlabFunction([(e1.'*e2 .* e3.' *e4),(e1.'*e3 .* e2.' *e4),...
    e1.'*e4 .* e2.' *e3]);
angle_fn = zeros(length(r_range),4,length(phi_range));
shift_fit = zeros(4,length(phi_range)); 
%figure; hold on
for lp1= 1:length(phi_range)
    phi = phi_range(lp1);
    % init_rot = [(sin(phi) + 1)/cos(phi),0,pi,-(cos(phi) + 1)/sin(phi)];  
     init_rot = [pi/2+phi,0,pi,pi+phi];  
shift = squeeze(minimum(:,lp1,:));%-repmat(minimum{1},length(b_range),1);
rnj = fval(lp1,:)< 2e-6 ; %good points;
%figure;plot(b_range(rnj),[shift(rnj,:),shift2(rnj,:)])

%plot(r_range(rnj),fval(lp1,rnj))
for lp2= 1:4
    
   tmp = ezfit( r_range(rnj).',shift(lp2,rnj),'a*x^2');%+b*x+c');
   shift_fit(lp2,lp1,:) =  tmp.m; 
   angle_fn(:,lp2,lp1) =  tmp.m(1).*r_range.^2;%+tmp.m(2).*r_range+tmp.m(3);
end

    for lp2 = 1:length(r_range)

        tmp = angle_fn(lp2,:,lp1)+init_rot;
        r_small = r_range(lp2)*cos(phi); %smaller one
    %min_smooth(lp) = fff(b_range(lp),minimum{1}(1)+shift_smooth(lp,1),minimum{1}(2)+shift_smooth(lp,2),...
    %minimum{1}(3)+shift_smooth(lp,3),minimum{1}(4)+shift_smooth(lp,4));

    min_smooth2(:,lp1,lp2) = fff(r_range(lp2),r_small,phi,tmp(1),tmp(2),tmp(3),tmp(4));  
    min_smooth(lp1,lp2) = f(r_range(lp2),r_small,phi,tmp(1),tmp(2),tmp(3),tmp(4));                      
                        
    end
end
% for lp = 1:length(b_long)
%     
%     min_smooth2(lp) = fff(b_long(lp),minimum{1}(1)+shift_long(lp,1),minimum{1}(2)+shift_long(lp,2),...
%     minimum{1}(3)+shift_long(lp,3),minimum{1}(4)+shift_long(lp,4));
% end

figure('Renderer','zbuffer'); pcolor(r_range,phi_range,min_smooth)
figure; plot(r_range,angle_fn(:,:,end))


%%

for lp = 2:length(b_range)
e1min = double(subs(ee1,[cc,b1],{minimum2(lp,1,3),b_range(lp)}));
e2min = double(subs(ee2,[cc,b1],{minimum2(lp,2,3),b_range(lp)}));
e3min = double(subs(ee3,[cc,b1],{minimum2(lp,3,3),b_range(lp)}));
e4min = double(subs(ee4,[cc,b1],{minimum2(lp,4,3),b_range(lp)}));
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