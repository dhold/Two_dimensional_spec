%finds angles of rotation related to faces of platonic solids

phi = (1+sym(sqrt(5)))/2;  %golden ratio
%list verticies of a regular icosahedron
vert_pos1 = [0,1,phi;0,-1,phi;0,1,-phi;0,-1,-phi]; %cyclic perms of these
vert_pos1 = [vert_pos1;vert_pos1(:,[3,1,2]);vert_pos1(:,[2,3,1])]/sqrt(1+phi^2); 

%list verticies of a regular dodecehedron   
vert_pos = sym([1,1,1;-1,1,1;1,-1,1;-1,-1,1]);
vert_pos = [vert_pos;-vert_pos];
v_tmp = [0,1/phi,phi; 0,-1/phi,phi ; 0,1/phi,-phi ; 0,-1/phi,-phi];
vert_pos = [vert_pos;v_tmp;v_tmp(:,[3,1,2]);v_tmp(:,[2,3,1])];
vert_pos = vert_pos/sqrt(3);
% 
% x =[1,0,0]; y=[0,1,0]; z =[0,0,1]; 

syms a b c real
%k vec along z direction
rz = [cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
ry = [cos(b),0,sin(b);0,1,0;-sin(b),0,cos(b)];
rz2 = [cos(c),-sin(c),0;sin(c),cos(c),0;0,0,1];
RR = rz2*ry*rz;
clear rot_set_iso rot_set_dodec

% rot_set_iso{1} = eye(3);
% for lp = 2:length(vert_pos1)
%     tmp = solve((vert_pos1(lp,:).')==RR*(vert_pos1(lp-1,:).'));
%     rot_set_iso{lp} = double(subs(RR,[a,b,c],[tmp.a,tmp.b,tmp.c]));
% end

% rot_set_dodec{1} = eye(3);
% for lp = 2:length(vert_pos)
%         tmp = solve((vert_pos(lp,:).')==RR*(vert_pos(lp-1,:).'));
%     rot_set_dodec{lp} = double(subs(RR,[a,b,c],[tmp.a,tmp.b,tmp.c]));
% end

%%
%{
   clear rot_set_iso_num  rot_set_dodec_num
rot_set_iso_num{1} = eye(3);  rot_set_dodec_num{1} = eye(3); 
    RRfun = @(a) [cos(a(3)),-sin(a(3)),0;sin(a(3)),cos(a(3)),0;0,0,1]*...
    [cos(a(2)),0,sin(a(2));0,1,0;-sin(a(2)),0,cos(a(2))]*...
        [cos(a(1)),-sin(a(1)),0;sin(a(1)),cos(a(1)),0;0,0,1];
    
   options.TolX = 1e-8;
   options.TolFun = 1e-8;

 for lp = 2:length(vert_pos1)
     
     v = double(vert_pos1(lp,:).'); vv =  double(vert_pos1(lp-1,:).');
    Rfun = @(a) norm(v-[cos(a(3)),-sin(a(3)),0;sin(a(3)),cos(a(3)),0;0,0,1]*...
    [cos(a(2)),0,sin(a(2));0,1,0;-sin(a(2)),0,cos(a(2))]*...
        [cos(a(1)),-sin(a(1)),0;sin(a(1)),cos(a(1)),0;0,0,1]*vv);     
     [A,B] = fminsearch(Rfun,[0,0,0],options);
     rot_set_iso_num{lp} = RRfun(A); 
 end
 
  for lp = 2:length(vert_pos)
     v = double(vert_pos(lp,:).'); vv =  double(vert_pos(lp-1,:).');
    Rfun = @(a) norm(v-[cos(a(3)),-sin(a(3)),0;sin(a(3)),cos(a(3)),0;0,0,1]*...
    [cos(a(2)),0,sin(a(2));0,1,0;-sin(a(2)),0,cos(a(2))]*...
        [cos(a(1)),-sin(a(1)),0;sin(a(1)),cos(a(1)),0;0,0,1]*vv);     
     [A,B] = fminsearch(Rfun,[0,0,0],options);
     rot_set_dodec_num{lp} = RRfun(A); 
 end
 save('Parameters/platonic_solid_rotations.mat','rot_set_iso_num','rot_set_dodec_num')
%}
%newer verison
   clear rot_set_iso_num  rot_set_dodec_num
rot_set_iso_num{1} = eye(3);  rot_set_dodec_num{1} = eye(3); 
    RRfun = @(a) [cos(a(3)),-sin(a(3)),0;sin(a(3)),cos(a(3)),0;0,0,1]*...
    [cos(a(2)),0,sin(a(2));0,1,0;-sin(a(2)),0,cos(a(2))]*...
        [cos(a(1)),-sin(a(1)),0;sin(a(1)),cos(a(1)),0;0,0,1];
    
   options.TolX = 1e-8;
   options.TolFun = 1e-8;
   vv = double(vert_pos1(1,:)).';

 for lp = 1:length(vert_pos1)
     
     v = double(vert_pos1(lp,:)).';
    Rfun = @(a) norm(v-[cos(a(3)),-sin(a(3)),0;sin(a(3)),cos(a(3)),0;0,0,1]*...
    [cos(a(2)),0,sin(a(2));0,1,0;-sin(a(2)),0,cos(a(2))]*...
        [cos(a(1)),-sin(a(1)),0;sin(a(1)),cos(a(1)),0;0,0,1]*vv);     
     [A,B] = fminsearch(Rfun,[0,0,0],options);
     rot_set_iso_num{lp} = RRfun(A); 
 end
 
 vv = double(vert_pos(1,:)).';
 
  for lp = 1:length(vert_pos)
     v = double(vert_pos(lp,:).'); 
    Rfun = @(a) norm(v-[cos(a(3)),-sin(a(3)),0;sin(a(3)),cos(a(3)),0;0,0,1]*...
    [cos(a(2)),0,sin(a(2));0,1,0;-sin(a(2)),0,cos(a(2))]*...
        [cos(a(1)),-sin(a(1)),0;sin(a(1)),cos(a(1)),0;0,0,1]*vv);     
     [A,B] = fminsearch(Rfun,[0,0,0],options);
     rot_set_dodec_num{lp} = RRfun(A); 
 end
 %save('Parameters/platonic_solid_rotations.mat','rot_set_iso_num','rot_set_dodec_num')
%%


rot_set_iso{1} = eye(3);  rot_set_dodec{1} = eye(3); 
 for lp = 2:length(vert_pos1)
     tmp = solve((vert_pos1(lp,:).')==RR*vert_pos1(lp-1,:).'); %will have 2 solns
     rot_set_iso{lp} = double(subs(RR,[a,b,c],[tmp.a(1),tmp.b(1),tmp.c(1)])); %pick first
 end

 for lp = 2:length(vert_pos)
         tmp = solve((vert_pos(lp,:).')==RR*vert_pos(lp-1,:).');
     rot_set_dodec{lp} = double(subs(RR,[a,b,c],[tmp.a(1),tmp.b(1),tmp.c(1)]));
 end

 save('Parameters/platonic_solid_rotations2.mat','rot_set_iso','rot_set_dodec')

%% DO it all symbolic innit



init_vec = [0;0;1]; %unit vector along z
R = rz*ry;
phi = (1+sym(sqrt(5)))/2;

vert_pos1 = [0,1,phi;0,-1,phi;0,1,-phi;0,-1,-phi]; %cyclic perms of these
vert_pos1 = [vert_pos1;vert_pos1(:,[3,1,2]);vert_pos1(:,[2,3,1])]/sqrt(1+phi^2); 

%list verticies of a regular dodecehedron   
vert_pos = [1,1,1;-1,1,1;1,-1,1;-1,-1,1];
vert_pos = [vert_pos;-vert_pos];
v_tmp = [0,1/phi,phi; 0,-1/phi,phi ; 0,1/phi,-phi ; 0,-1/phi,-phi];
vert_pos = [vert_pos;v_tmp;v_tmp(:,[3,1,2]);v_tmp(:,[2,3,1])];
vert_pos = vert_pos/sqrt(3);


 for lp = 1:length(vert_pos1)
     tmp = solve((vert_pos1(lp,:).')==R*init_vec); %will have 2 solns
     rot_set_iso{lp} = double(subs(R,[a,b],[tmp.a(1),tmp.b(1)])); %pick first
 end

 for lp = 1:length(vert_pos)
         tmp = solve((vert_pos(lp,:).')==R*init_vec);
     rot_set_dodec{lp} = double(subs(R,[a,b],[tmp.a(1),tmp.b(1)]));
 end



 %%
 save('Parameters/platonic_solid_rotations3.mat','rot_set_iso','rot_set_dodec')