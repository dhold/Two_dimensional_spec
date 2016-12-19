
N = 3; H_e = [1000,100,50;100,1200,150;50,150,1300];
%N = 2; H_e = [1000,100;100,1200];
[H_el,H_vib,H_e,H_f,e_1,e_2,fock_space_rep,M] = ...
    generate_ex_vib_ham3(H_e,[],[],[],[],[],[],[]) ;

mu = rand(N,3); R = rand(N,3); mm =zeros(N,3);

[C,~] = eig(H_el(2:N+1,2:N+1));
[CC,~] = eig(H_el(N+2:end,N+2:end));
fock2 = fock_space_rep(N+2:end,:);

mu_e_f = zeros(N*(N-1)/2,3,N); mu_ex = zeros(N,3);
for k = 1:N
    for n = 1:N
    mu_ex(k,:) = mu_ex(k,:) + C(n,k)*mu(n,:); %this is what I think is
    %correct as I need the participation of the nth site in the kth exciton
    %the wavefunction for the kth exciton is 
    %|psi_k> = sum_n C(n,k) |psi_n>
    %mu_ex(k,:) = mu_ex(k,:) + C(k,n)*mu(n,:);
    end
    for f = 1:N*(N-1)/2

for n=1:N
    for m = n+1:N
        lg = logical(fock2(:,n))&logical(fock2(:,m));
        mu_e_f(f,:,k) =  mu_e_f(f,:,k) + CC(lg,f)*(C(n,k)*mu(m,:) + C(m,k)*mu(n,:));
       % mu_e_f(f,:,k) =  mu_e_f(f,:,k) + CC(f,lg)*(C(k,n)*mu(m,:) + C(k,m)*mu(n,:)); 
    end
end
    end
end

MM = blkdiag(1,M{1},M{2});

V_n = zeros(1+N+N*(N-1)/2,1+N+N*(N-1)/2,N); sz2 = length(H_vib);

for lp =1:N
   
    V = zeros(1+N+N*(N-1)/2);  
    V(1,lp+1) = 1;  %mixes to ground
    
    lg = [false(N+1,1);fock_space_rep(N+2:end,lp)==1]; %elements with excitation at lp
    lg2 = fock_space_rep; lg2(:,lp) = 0; %[~,lg2] = find(lg2(lg,:)); 
    V(lg,2:N+1) =lg2(lg,:);
    V = V+V'; V_n(:,:,lp) = V;
    
end


param_set = [[0,0,1];[0,0,1];[0,1,0];[1,0,0];[1,0,0];[1,0,0];[1,0,0]];
[~,~,~,av_4,av_5,av2_5]=ori_precomp_site_basis(param_set,mu,R,imag(mm));
%also compute some first order averages with CP light
 tmp = [0,0,1;[1,+1i,0]/sqrt(2);[1,-1i,0]/sqrt(2)];
 
[av_2,av_3,av2_3]=ori_precomp_site_basis(tmp,mu,R,imag(mm));
%project exciton basis interaction operator into the site basis
alpha_n = mtimesx(MM,mtimesx(V_n,MM,'C')); 
Cg  = squeeze(alpha_n(1,2:N+1,:)) ;  Cf = alpha_n(2:N+1,N+2:end,:);
%%

%1st order averages
av_set_fo = mtimesx(mtimesx(Cg,'C',av_2),Cg); 
av_set_cd = mtimesx(mtimesx(Cg,'C',av_3+av2_3),Cg); 

%3rd order averages

tmp = mtimesx(mtimesx(Cg,'C',av_4),Cg); %first two transitions are between
%g and one exciton transitions

%the next loop is more complicated because it involves transitions from e-f
%and then f'-e', need to reshape
av_set_ESA = zeros(N,N,N,N,N*(N-1)/2,N*(N-1)/2,size(av_set_GSB,5));
sz = size(tmp); 
if ndims(tmp) ==4
    sz = [sz,1];
end

%Terrible giant for loop
tic
for k3=1:N
    for k4=1:N
for f1 = 1:N*(N-1)/2
    for f2 = 1:N*(N-1)/2
        tmp2 = zeros(sz(1),sz(1),1,1,sz(5));
        for flp = 1:N*(N-1)/2
            ss = find(fock2(flp,:));
            for flp2 = 1:N*(N-1)/2
                ss2 = find(fock2(flp2,:));
                
tmp2 = tmp2 + M{2}(flp,f1)*M{2}(flp2,f2)*(...
        M{1}(ss(1),k3)*M{1}(ss2(1),k4)*tmp(:,:,ss(2),ss2(2),:) +...
        M{1}(ss(2),k3)*M{1}(ss2(1),k4)*tmp(:,:,ss(1),ss2(2),:) +...
        M{1}(ss(1),k3)*M{1}(ss2(2),k4)*tmp(:,:,ss(2),ss2(1),:) +...
        M{1}(ss(2),k3)*M{1}(ss2(2),k4)*tmp(:,:,ss(1),ss2(1),:));

            end
        end
        av_set_ESA(:,:,k3,k4,f1,f2,:) = tmp2;        
    end
end
    end
end
toc

tmp = permute(tmp,[3,4,1,2,5]); %permute so next two dimensions are in mtimesx
av_set_ESA_test =  zeros(N,N,N,N,N*(N-1)/2,N*(N-1)/2,size(av_set_GSB,5));
tic
for k3=1:N %much bette loop
    for k4=1:N
        tmp2 = zeros(N*(N-1)/2,N*(N-1)/2,sz(1),sz(2),sz(5));
        for flp = 1:N*(N-1)/2
            ss = find(fock2(flp,:));
            for flp2 = 1:N*(N-1)/2
                ss2 = find(fock2(flp2,:));
                
tmp2 = tmp2 + mtimesx(M{2}(flp,:).'*M{2}(flp2,:),(...
        M{1}(ss(1),k3)*M{1}(ss2(1),k4)*tmp(ss(2),ss2(2),:,:,:) +...
        M{1}(ss(2),k3)*M{1}(ss2(1),k4)*tmp(ss(1),ss2(2),:,:,:) +...
        M{1}(ss(1),k3)*M{1}(ss2(2),k4)*tmp(ss(2),ss2(1),:,:,:) +...
        M{1}(ss(2),k3)*M{1}(ss2(2),k4)*tmp(ss(1),ss2(1),:,:,:)));

            end
        end
        av_set_ESA_test(:,:,k3,k4,:,:,:) = permute(tmp2,[3,4,1,2,5]);        
    end
end
toc


tmp = squeeze(tmp);
%the first one involves transitions between g-e only (GSB, SE)
av_set_GSB = mtimesx(mtimesx(Cg,'C',tmp),Cg);
av_set_GSB = permute(av_set_GSB,[3,4,1,2,5]); %permute back to original order

% Cf2=permute(Cf,[3,1,2]); %reshape this way
% 
% av_set_ESA = zeros(N,N,N,N,N*(N-1)/2,N*(N-1)/2,size(av_set_GSB,5));
% for f = 1:N*(N-1)/2
%     for f2 = 1:N*(N-1)/2
% tmp2 = mtimesx(mtimesx(Cf2(:,:,f),'C',tmp),Cf2(:,:,f2)); %ESA type
% tmp2 = permute(tmp2,[3,4,1,2,5]);
% av_set_ESA(:,:,:,:,f,f2,:) = tmp2;
%     end
% end
%av_set_ESA (e1,e2,e3*(N-1)*N/2+f1,e4*(N-1)*N/2+f2) type calls
%next calculate the higher order moments from the other interactions

tmp = mtimesx(mtimesx(Cg,'C',av_5+av2_5),Cg);
tmp = permute(tmp,[3,4,1,2,5,6]); 


%%
pol_linear = [[1,0,0];[1,0,0];[1,0,0];[1,0,0]];
av_lin_pol = zeros(N,N,N,N);
av_lin_pol_f = zeros(N,N,N,N,N*(N-1)/2,N*(N-1)/2);
for k1 = 1:N %slow loop
    for k2 = 1:N
       for k3 = 1:N 
           for k4 = 1:N  
              
              mu_set = [mu_ex(k1,:);mu_ex(k2,:);mu_ex(k3,:);mu_ex(k4,:)];
              av_lin_pol(k1,k2,k3,k4) = tensor_av(mu_set,pol_linear);
              for f=1:N*(N-1)/2
                  for f2=1:N*(N-1)/2
                  mu_set = [mu_ex(k1,:);mu_ex(k2,:);mu_e_f(f,:,k3);mu_e_f(f2,:,k4)];
                  av_lin_pol_f(k1,k2,k3,k4,f,f2) = tensor_av(mu_set,pol_linear);
                  end
              end
           end
       end
    end
end

%% 
test1 = av_lin_pol-av_set_GSB; sum(abs(test1(:)))
test2 = av_lin_pol_f-av_set_ESA; sum(abs(test2(:)))
test3 = av_set_ESA-av_set_ESA_test;  sum(abs(test3(:)))