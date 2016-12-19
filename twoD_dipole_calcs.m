%% Solve to find vectors satisfying the condition
syms theta1 theta2 theta3 real
%theta4 = 0
[sol1,sol2,sol3] = solve(0==cos(theta3)*cos(theta2-theta1),...
-1/2==cos(theta2)*cos(theta3-theta1),...
1/2==cos(theta1)*cos(theta3-theta2));
%% check
j=1
double(0-cos(sol3(j))*cos(sol2(j)-sol1(j)))
double(-1/2-cos(sol2(j))*cos(sol3(j)-sol1(j)))
double(1/2-cos(sol1(j))*cos(sol3(j)-sol2(j)))


%%
M4 = 1/30*[4,-1,-1;-1,4,-1;-1,-1,4];
syms T1 T2 T3 t1 t2 t3 t4 TT1 TT2 TT3 TT4 tt1 tt2 tt3 tt4 real
% T4 = zero (angles relative to this beam) ,
% t -> dipole angles; TT -> magfield angles; tt -> magdipole angles
ori_av = [cos(TT4-T3)*cos(T2-T1),cos(TT4-T2)*cos(T3-T1),cos(TT4-T1)*cos(T3-T2)]...
   *M4*[cos(tt4-t3)*cos(t2-t1);cos(tt4-t2)*cos(t3-t1);cos(tt4-t1)*cos(t3-t2)]  ;

%% General magic angle configuation

syms theta1 theta2 theta3 C
%%
[sol1,sol2] = solve(...
cos(theta2)*cos(theta1)==cos(theta2-theta1)/3,...
cos(theta1)*cos(theta2)==cos(theta2-theta1)/3);
%%
[sol1,sol2,sol3] = solve(...
cos(theta2)*cos(theta1-theta3)==cos(theta3)*cos(theta2-theta1)/3,...
cos(theta1)*cos(theta2-theta3)==cos(theta3)*cos(theta2-theta1)/3);
%%
[sol11,sol22,sol33,solC] = solve(...
cos(theta2)*cos(theta1-theta3)==cos(theta3)*cos(theta2-theta1)/3+6*C,...
cos(theta1)*cos(theta2-theta3)==cos(theta3)*cos(theta2-theta1)/3-6*C);

%%
M5 = 1/30*[3,-1,-1,1,1,0;-1,3,-1,-1,0,1;-1,-1,3,0,-1,-1;...
        1,-1,0,3,-1,1;1,0,-1,-1,3,-1;0,1,-1,1,-1,3];
    syms A B C D E F 
    some_vec = [A,B,C,D,E,F];
    
sol = solve( some_vec*M5==[1,0,0,0,0,0]);   
ans1 = [sol.A,sol.B,sol.C,sol.D,sol.E,sol.F];
sol = solve( some_vec*M5==[0,1,0,0,0,0]);   
ans2 = [sol.A,sol.B,sol.C,sol.D,sol.E,sol.F];
sol = solve( some_vec*M5==[0,0,1,0,0,0]);   
ans3 = [sol.A,sol.B,sol.C,sol.D,sol.E,sol.F];
sol = solve( some_vec*M5==[0,0,0,1,0,0]);   
ans4 = [sol.A,sol.B,sol.C,sol.D,sol.E,sol.F];
sol = solve( some_vec*M5==[0,0,0,0,1,0]);   
ans5 = [sol.A,sol.B,sol.C,sol.D,sol.E,sol.F];
sol = solve( some_vec*M5==[0,0,0,0,0,1]);   
ans6 = [sol.A,sol.B,sol.C,sol.D,sol.E,sol.F];

ans_tot = [ans1;ans2;ans3;ans4;ans5;ans6]


%% Calculate the factors for a general setup in terms of angles

%uses that A_1 x A_2 . A_3 = pm
%sqrt(S_12^2 S_13^2 +2 C_13 C_23 C_12 - C_23^2- C_12^2 C_13 ^2)
%for real unit vectors A_1 A_2 A_3 and C_ij = Cos(theta_ij), S_ij =
%sin(theta_ij) etc.
M4 = 1/30*[4,-1,-1;-1,4,-1;-1,-1,4];
M5 = 1/30*[3,-1,-1,1,1,0;-1,3,-1,-1,0,1;-1,-1,3,0,-1,-1;...
        1,-1,0,3,-1,1;1,0,-1,-1,3,-1;0,1,-1,1,-1,3];
    
%syms t_12 t_13 t_14 t_15 t_23 t_24 t_25 t_35 t_34 t_45 real

F5_sym = @(t_12, t_13, t_14, t_15, t_23, t_24, t_25, t_34, t_35, t_45) 1/30 * ...
   [cos(t_15)*sqrt(sin(t_23)^2*sin(t_34)^2 - (cos(t_24)-cos(t_34)*cos(t_23))^2),...
    cos(t_25)*sqrt(sin(t_13)^2*sin(t_34)^2 - (cos(t_14)-cos(t_34)*cos(t_13))^2),... %1->2 from 1
    cos(t_12)*sqrt(sin(t_35)^2*sin(t_34)^2 - (cos(t_45)-cos(t_34)*cos(t_35))^2),... %5->1 from 2
    cos(t_35)*sqrt(sin(t_12)^2*sin(t_24)^2 - (cos(t_14)-cos(t_24)*cos(t_12))^2),... %2->3 from 2
    cos(t_13)*sqrt(sin(t_25)^2*sin(t_24)^2 - (cos(t_45)-cos(t_24)*cos(t_25))^2),... %2->3 from 3
    cos(t_23)*sqrt(sin(t_15)^2*sin(t_14)^2 - (cos(t_45)-cos(t_14)*cos(t_15))^2)];   %1->2 from 5


