function [R1,R2,R3,R4,R5,R6] = rp_cont_HEOM2(sz_params,rho_eg,V_ge_t,V_ef_t,...
                            prop_op_gg,prop_op_ee,t2_range,e1,e2,...
      calc_rp,calc_nr,av1_g,av1_f,av2_g,av2_f,H_trunc,tol) ;

%Calculated the re/nonrephasing contribution to a 2D signal given the coherence
%part of the density matrix, the Heisenberg picture operators V(t_3) and 
% the operators that propogate the system in the ground and excited state
% This will sum over all pathways with the dipole
% averaging constants in av1, av2 to give a result.  It is assumed that the
% av2 are of the averaging form in my notes with the final index running
% from 1-3 each of which is weighted by |k1| |k2| |k3| etc and the factor
% of +/- 1 depending on which signal we are considering

%sz_params is a vector with the important length constants, note that code
%works (should) even if HL = 1 (not the HEOM, just a normal propogator)
N = sz_params(1); %number of sites
HL = sz_params(2); %total number of aux density matricies (+1)
t1L = size(rho_eg,4);  t2L = length(t2_range);  t3L = size(V_ge_t,2);
s_rp = [-1,+1,+1]; s_nr = [+1,-1,+1];

if ~exist('tol','var')
    tol =[1e-7,1e-5];
end
out_trunc_gg = 1*H_trunc;
out_trunc_ee = N^2*1*H_trunc;


if ~calc_rp && ~calc_nr
   warning('nothing calculated')
   R1=[]; R2=[]; R3=[]; R4=[]; R5=[]; R6=[]; return
end
    num_sets = [size(av1_g,5),size(av2_g,5)];
    %number of sets of values to calculate for the 4th order type averages
    %and 5th order type averages
    
    if calc_rp
 %preallocate rephasing
        R1 = zeros(t3L,t2L,t1L,sum(num_sets)); R4 = R1; R5 = R4;
    else
        R1=[]; R4=[]; R5=[];
    end
    if calc_nr
    R2 = zeros(t3L,t2L,t1L,sum(num_sets));  %preallocate nonrephasing
        R3 = R2; R6 = R3;
    else
        R2 = []; R3=[]; R6 = [];
    end    

rho_eg = reshape(rho_eg,N,HL,t1L,N); %reshape to this form for convinience

%reshape V_ex_eg and V_ex_ef to the appropriate sizes to act on the HEOM
%and put them in Lioville space
% for j = 1:N
%    %note this is the opposite direction to the over arrows, it is the side
%    %the operator acts from
%     V_eg_left{j} = repmat(sparse(  kron(V_ex_eg(:,:,j),eye(size(V_ex_eg(:,:,j))))),[HL,1]); %#ok<*AGROW>
%     V_eg_right{j} = repmat(sparse(  kron(eye(size(V_ex_eg(:,:,j))),V_ex_eg(:,:,j))),[HL,1]);
%     V_ge_left{j} = repmat(sparse(  kron(V_ex_eg(:,:,j)',eye(size(V_ex_eg(:,:,j)')))),[HL,1]);
%     V_ge_right{j} = repmat(sparse(  kron(eye(size(V_ex_eg(:,:,j)')),V_ex_eg(:,:,j)')),[HL,1]);
%     for f=1:N*(N-1)/2
%     V_fe_left{j,f} = repmat(sparse(  kron(V_ex_fe(:,:,j,f),eye(size(V_ex_fe(:,:,j,f))))),[HL,1]);    
%     end
% end

%generate maps which can pick out the hermitian conjugate of a Liouville
%space vector

temp = reshape(N^2,N,N); temp2 = temp';  cnj_ee = reshape(temp2,N^2,1);

 cj_ee = HEOM_conj({HL,N*1,N*1});
 cj_gg = HEOM_conj({HL,1,1});
 cj_eg = HEOM_conj({HL,N*1,1});
 cj_fe = HEOM_conj({HL,N*(N-1)/2*1,N*1});

rho_gg_prop = @(t,v) prop_op_gg*v;
rho_ee_prop = @(t,v) prop_op_ee*v;

for j1 = 1:N  %first transition to exciton j and then either conjugate or j->f
    for j2 = 1:N %2nd dipole interaction
        
        %propogate the resulting excited state pop and ground state hole
        
    tmp_gg = rho_eg(j2,:,:,j1);   %picks these elements
    tmp_gg = reshape(tmp_gg,[HL,t1L]); %reshape to normal Liouville space
    
    %reshape rho_eg into seperate Heirarchy components 
    tmp_ee =  repmat(rho_eg(:,:,:,j1),[N,1,1]);
    lg = true(N^2,1); lg(1+(j2-1)*N:j2*N) = false; %these elements remain
    tmp_ee(lg,:,:) = 0;

    tic
    %loop over t1L, yes this WILL be painfully slow if this is long
   % figure 
   % hold on
    for lp = 1:t1L %could avoid nexting other loops within this
        %but it might mean quite a lot of saved data.
        tmp_gg_sec = tmp_gg(:,lp);
        tmp_gg_t = (OD_wrapper(t2_range,rho_gg_prop,tmp_gg_sec,out_trunc_gg,'ode45',tol)).';      
        %tmp_gg_t = reshape(tmp_gg_t,[1,HL,t2L]); kind of pointless
        tmp_ee_sec = tmp_ee(:,lp);
        tmp_ee_t = (OD_wrapper(t2_range,rho_ee_prop,tmp_ee_sec,out_trunc_ee,'ode45',tol)).';
        tmp_ee_t = reshape(tmp_ee_t,[N^2,HL,t2L]); %reshape
       % plot(t2_range,tmp_ee_t([1,4],:))
        %now calculate each of the contributions
        for j3 = 1:N
            
            if calc_rp
           % R1tmp = %stimulated emission
           % R4tmp = %-as abs but negative occ
%            lg = false(N^2,1); lg(1+N*(j3-1):j3*N) = true; 
%            lg = repmat(lg,[HL,1]);
           
          R1tmp = tmp_ee_t(1+N*(j3-1):j3*N ,:,:);
          R4tmp = zeros(N*HL,t2L); R4tmp(j3:N:N*HL,:) = tmp_gg_t;
            end
             if calc_nr
           % R2tmp = V_eg_right{j3}*conj(tmp_ee_t(cj_ee,:)); %SE
           % R3tmp = -V_eg_left{j3}*conj(tmp_gg_t(cj_gg,:)); %-as abs
%            lg = false(N^2,1); lg(1+N*(j3-1):j3*N) = true; 
%            lg = repmat(lg,[HL,1]);
           
          R2tmp = conj(tmp_ee_t(cnj_ee,:,:)); %conjugate the matrix
          R2tmp = R2tmp(1+N*(j3-1):j3*N ,:); %take correct elements out of it
          %rho_gg is only size one and should be real so no transpose
          %is needed, but underdamped mode aux mat are complex so need conj
          R3tmp = zeros(N*HL,t2L); R3tmp(j3:N:N*HL,:) = conj(tmp_gg_t);
             end
            
            for j4 = 1:N       
              V_ge_sec = V_ge_t(:,:,j4);
            if calc_rp                  
                %first pure dipole factor
            ori_tmp1 = reshape(av1_g(j1,j2,j3,j4,:),[1,1,num_sets(1)]);
            %next the one which depends on |k_j|
            tmp2 = squeeze(av2_g(j1,j2,j3,j4,:,:));
            ori_tmp2 = reshape(e1(j1)*s_rp(1).*tmp2(:,1)+e1(j2)*s_rp(2).*tmp2(:,2)+...
                        e1(j3)*s_rp(3).*tmp2(:,3),[1,1,num_sets(2)]);
            ori_fct_rf = cat(3,ori_tmp1,ori_tmp2);
             %combine with the contributions, flip round t2 and t3    
            %this operation is essentially a trace

            R_1 = mtimesx( V_ge_sec,'C',R1tmp); %t3_range by t2_range
            R_4 = mtimesx( V_ge_sec,'C',R4tmp);

            R_1 = repmat(R_1,[1,1,length(ori_fct_rf)]).*repmat(ori_fct_rf,[t3L,t2L,1]);
            R_4 = repmat(R_4,[1,1,length(ori_fct_rf)]).*repmat(ori_fct_rf,[t3L,t2L,1]);
            end
            if calc_nr
            ori_tmp1 = reshape(av1_g(j1,j2,j3,j4,:),[1,1,num_sets(1)]);
            %next the one which depends on |k_j|
            tmp2 = squeeze(av2_g(j1,j2,j3,j4,:,:));
            ori_tmp2 = reshape(e1(j1)*s_nr(1).*tmp2(:,1)+e1(j2)*s_nr(2).*tmp2(:,2)+...
                        e1(j3)*s_nr(3).*tmp2(:,3),[1,1,num_sets(2)]);
            ori_fct_nr = cat(3,ori_tmp1,ori_tmp2);                

            R_2 =  mtimesx(V_ge_sec,'C',R2tmp); %t3_range by t2_range
            R_3 =  mtimesx(V_ge_sec,'C',R3tmp);
            R_2 = repmat(R_2,[1,1,length(ori_fct_rf)]).*repmat(ori_fct_nr,[t3L,t2L,1]);
            R_3 = repmat(R_3,[1,1,length(ori_fct_rf)]).*repmat(ori_fct_nr,[t3L,t2L,1]);
            end       
            if calc_rp
                R_6 = zeros(size(R_1));
            end
            if calc_nr
                R_5 = zeros(size(R_2));
            end
                        
                for f = 1:N*(N-1)/2
                    if calc_rp
                  % R6tmp = -(V_fe_left{j3,f}*tmp_ee_t); 

                 lg = true(N,N*(N-1)/2); lg(:,f) = false;
                 lg=reshape(lg,[N^2*(N-1)/2,1]); lg = repmat(lg,[HL,1]);
                 lg2 = false(N^2,1); lg2(1+N*(j3-1):N*j3)=true;
                 tmp_ee_sec = tmp_ee_t(lg2,:);
                 R6tmp = zeros(N^2*(N-1)/2*HL,t2L); 
                 R6tmp = zeros(N^2*(N-1)/2*HL,t2L);
                   R6tmp = -HEOM_mult(V_ex_fe(:,:,j3,f),'n',tmp_ee_t,'n','L',...
                        [N*1,N*1,HL]); 
             ori_tmp1 = reshape(av1_f(j1,j2,N*(N-1)/2*(j3-1)+f,N*(N-1)/2*(j4-1)+f,:),[1,1,num_sets(1)]);
            %next the one which depends on |k_j|
            tmp2 = squeeze(av2_f(j1,j2,N*(N-1)/2*(j3-1)+f,N*(N-1)/2*(j4-1)+f,:,:));
            ori_tmp2 = reshape(e1(j1)*s_rp(1).*tmp2(:,1)+e1(j2)*s_rp(2).*tmp2(:,2)+...
                        (e2(f)-e1(j3))*s_rp(3).*tmp2(:,3),[1,1,num_sets(2)]);
            ori_fct_rp = cat(3,ori_tmp1,ori_tmp2);                      
                   
            R6tmp2 = mtimesx( V_ef_t(cj_fe,:,j4,f),'T',R6tmp);
          % R6tmp2 = mtimesx( V_ef_t(:,:,j4,f),'T',R6tmp);
            R_6 = R_6 + repmat(R6tmp2,[1,1,length(ori_fct_rf)]).*repmat(ori_fct_rp,[t3L,t2L,1]);
              %combine with the contributions 
                    end
                    if calc_nr
                    %R5tmp = -V_fe_left{j3,f}*conj(tmp_ee_t(cj_ee,:));    
                    R5tmp = -HEOM_mult(V_ex_fe(:,:,j3,f),'n',tmp_ee_t,'C','L',...
                        [N*1,N*1,HL]); 
             ori_tmp1 = reshape(av1_f(j1,j2,N*(N-1)/2*(j3-1)+f,N*(N-1)/2*(j4-1)+f,:),[1,1,num_sets(1)]);
            %next the one which depends on |k_j|
            tmp2 = squeeze(av2_f(j1,j2,N*(N-1)/2*(j3-1)+f,N*(N-1)/2*(j4-1)+f,:,:));
            ori_tmp2 = reshape(e1(j1)*s_nr(1).*tmp2(:,1)+e1(j2)*s_nr(2).*tmp2(:,2)+...
                        (e2(f)-e1(j3))*s_nr(3).*tmp2(:,3),[1,1,num_sets(2)]);
            ori_fct_nr = cat(3,ori_tmp1,ori_tmp2);       
            
           R5tmp2 = mtimesx( V_ef_t(cj_fe,:,j4,f),'T',R5tmp);
           %R5tmp2 = mtimesx( V_ef_t(:,:,j4,f),'T',R5tmp);
            R_5 = R_5 + repmat(R5tmp2,[1,1,length(ori_fct_rf)]).*repmat(ori_fct_nr,[t3L,t2L,1]); 
                   
                    end                    
                end
   %combine this frequency slices into the total response function             
               if calc_rp   
                   for lp2 = 1:size(R_2,3)
      R2(:,:,lp,lp2) = R_2(:,:,lp2); R3(:,:,lp,lp2) = R_3(:,:,lp2); 
      R5(:,:,lp,lp2) = R_5(:,:,lp2);
                   end
               end
               if calc_nr
                  for lp2 = 1:size(R_1,3)
      R1(:,:,lp,lp2) = R_1(:,:,lp2); R4(:,:,lp,lp2) = R_4(:,:,lp2); 
      R6(:,:,lp,lp2) = R_6(:,:,lp2);
                  end
               end
            end
        end
    end
    toc
    end
end
end



function rho_out = HEOM_mult(V,m1,rho_vec,m2,m3,sz)
%multiplies each tier of the HEOM matrix set by the vector V, arguments
%m1,m2 must be strings, 'n' for normal 'T' for transpose etc, see mtimesx
%m3 determines what the order of multiplication is. sz is a length two
%vector which determines what size blocks 

sz_act = size(rho_vec);

rho_out = reshape(rho_vec,[sz,sz_act(2)]);

if m3 == 'L'
rho_out = mtimesx(V,m1,rho_out,m2);
elseif  m3 == 'R'
rho_out = mtimesx(rho_out,m2,V,m1);   
end

rho_out = reshape(rho_out,[numel(rho_out)/sz_act(2),sz_act(2)]);


end

function rho_conj = HEOM_conj(rho_vec)
%conjugate transposes each block in the HEOM 

persistent conj_map

if iscell(rho_vec) %generate map which transposes
   
    HL_tmp = rho_vec{1}; N1_tmp = rho_vec{2}; N2_tmp = rho_vec{3};
    %create the logical map which pics the elements which will be
    %conjugated
    
    temp = reshape(1:N1_tmp*N2_tmp,N1_tmp,N2_tmp); temp2 = temp';
    temp3 = reshape(temp2,N1_tmp*N2_tmp,1); %new mappings
    
    conj_map = zeros(N1_tmp*N2_tmp*HL_tmp,1);
    for jj = 1:HL_tmp
        
        conj_map((jj-1)*N1_tmp*N2_tmp+1:jj*N1_tmp*N2_tmp) = temp3 + (jj-1)*N1_tmp*N2_tmp;
        
    end
    rho_conj = conj_map; return
end
    
    rho_conj = conj(rho_vec(conj_map));


end