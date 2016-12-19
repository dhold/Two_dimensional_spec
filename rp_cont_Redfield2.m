function [R1,R2,R3,R4,R5,R6,rho_gg_t,rho_ee_t] = rp_cont_Redfield2(sz_params,...
            rho_eg,rho_ge,V_ge,V_ef,prop_op_gg,prop_op_ee,t2_range,e1,e2,...
      calc_nr,calc_rp,av1_g,av1_f,av2_g,av2_f,tol,t_verbose) ;

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
N2 = N*(N-1)/2; %number of double excited states\
Nv = sz_params(2); %number of vibrational states included
%lengths of the three different times
t1L = size(rho_eg,2);  t2L = length(t2_range);  t3L = size(V_ge,1);

%order of conjugation of pulses for each signal, 
%this is the sign in e^(+/- i(k dot r-omega t)), in order they interact
s_rp = [-1,+1,+1];%[+1,-1,+1]; %rephasing or photon echo
s_nr = [+1,-1,+1]; %non rephasing signal

pause on %nice to be able to stop the code in windows

if ~exist('tol','var')
    tol =[1e-7,1e-5];
elseif isempty(tol)
    tol =[1e-7,1e-5];
end
if ~exist('t_verbose','var')
    t_verbose = false; %don't output time for loops
end

if nargout > 6
   rho_gg_t = zeros(Nv^2,t2L,t1L,N,N);
   rho_ee_t  = zeros(N^2*Nv^2,t2L,t1L,N,N);
end
out_trunc_gg = [];
out_trunc_ee = [];

if ~calc_nr && ~calc_rp
   warning('nothing calculated')
   R1=[]; R2=[]; R3=[]; R4=[]; R5=[]; R6=[]; return
end
    num_sets = [size(av1_g,5),size(av2_g,5)];
    %number of sets of values to calculate for the 4th order type averages
    %and 5th order type averages
    
    if calc_nr
 %preallocate rephasing
        R1 = complex(zeros(t3L,t2L,t1L,sum(num_sets))); R4 = R1; R5 = R4;
    else
        R1=[]; R4=[]; R5=[];
    end
    if calc_rp
    R2 = complex(zeros(t3L,t2L,t1L,sum(num_sets)));  %preallocate nonrephasing
        R3 = R2; R6 = R3;
    else
        R2 = []; R3=[]; R6 = [];
    end    

%reshape V_eg and V_ef to the appropriate sizes to act on the HEOM
%and put them in Lioville space.  This is different when they act on
%different components
%1) reshape(A * C, NN_1*NN_2,1) = kron(eye(NN_2),A)*reshape(C, NN_1*NN_2,1)
%2) reshape(C * A, NN_1*NN_2,1) = kron(A.',eye(NN_1))*reshape(C,NN_1*NN_2,1)
% For N_2 by N_1 A and N_2 by N_3 C we have
% reshape(A * C, N_1*N_3,1) == kron(eye(N_3),A)*reshape(C,N_3*N_2,1)
% For N_1 by N_2 A and N_3 by N_2 C we have
% reshape(C * A, N_1*N_3,1) == kron(A.',eye(N_3))*reshape(C,N_3*N_2,1)
for j = 1:N
   %this uses the arrow convention of "Advancing HEOM for efficient
   %evalution of ..." 1) is a right arrow 2) is a left arrow
   %V_eg = zeros(N,1);  V_ge = zeros(1,N); V_eg(j) = 1;  V_ge(j) = 1; 
    VV_eg = sparse((j-1)*Nv+1:j*Nv,1:Nv,ones(Nv,1),N*Nv,Nv); 
    VV_ge = sparse(1:Nv,(j-1)*Nv+1:j*Nv,ones(Nv,1),Nv,N*Nv);

   V_eg_gg_right{j} = kron(speye(Nv),VV_eg);     %#ok<*AGROW>
    
   V_eg_right{j} = kron(speye(N*Nv),VV_eg);
   V_eg_left{j} = kron(VV_eg.',speye(Nv));     
   V_eg_ee_left{j} = kron(VV_eg.',speye(N*Nv));
   
   V_ge_right{j} = kron(speye(Nv),VV_ge);
   V_ge_left{j} = kron(VV_ge.',speye(N*Nv));        

    for f3=1:N2
        VV_fe = sparse((f3-1)*Nv+1:f3*Nv,(j-1)*Nv+1:j*Nv,ones(Nv,1),N2*Nv,N*Nv);
        V_fe_right{j,f3} = kron(speye(N*Nv),VV_fe);
    end
end

rho_gg_prop = @(t,v) prop_op_gg*v;
rho_ee_prop = @(t,v) prop_op_ee*v;

for j1 = 1:N  %first transition to exciton j and then either conjugate or j->f
    %for j2 = j1 %2nd dipole interaction
        %in two colour experiment with om_1 = om_2 others are heavily
        %suppressed to the point they can be neglected
     for j2 = 1:N %2nd dipole interaction   
        %Calculate resulting excited state pop and ground state hole
        
   % NR-> rho_ee = rho__eg V_ge, rho_gg = V_ge rho__eg        
    tmp_gg = V_ge_right{j2}*rho_eg(:,:,j1); 
    tmp_ee =  V_ge_left{j2}*rho_eg(:,:,j1); 
         %might be best to use some logic here that if the element is too
         %small to bother with I just set the range for all t2/t3 to zero
         %lg = trace(tmp_gg)>thresh;
         %lg2 = trace(tmp_ee)>thresh;
         
         
    % RP-> rho_ee = V_eg rho__ge , rho_gg = rho__ge V_eg 
    tmp_gg2 = V_eg_left{j2}*rho_ge(:,:,j1);   %from conjugate path
    tmp_ee2 = V_eg_right{j2}*rho_ge(:,:,j1); 
if t_verbose
    tic
end
    %loop over t1L, yes this WILL be painfully slow if t1L is large and 
    % t2_range is long.  This is essentially unavoidable
pause(1e-3) %take a millisecond pause to allow for interupting the code
    for lp = 1:t1L %could avoid nexting other loops within this
        %but it might mean quite a lot of saved data.

        if calc_nr
        tmp_gg_sec = tmp_gg(:,lp);
        tmp_gg_t = (OD_wrapper(t2_range,rho_gg_prop,tmp_gg_sec,out_trunc_gg,'ode45',tol)).';      
        tmp_ee_sec = tmp_ee(:,lp);
        tmp_ee_t = (OD_wrapper(t2_range,rho_ee_prop,tmp_ee_sec,out_trunc_ee,'ode45',tol)).';
        end
        if calc_rp
        tmp_gg_sec2 = tmp_gg2(:,lp);
        tmp_gg_cnj = (OD_wrapper(t2_range,rho_gg_prop,tmp_gg_sec2,out_trunc_gg,'ode45',tol)).';      

        tmp_ee_sec2 = tmp_ee2(:,lp);
        tmp_ee_cnj = (OD_wrapper(t2_range,rho_ee_prop,tmp_ee_sec2,out_trunc_ee,'ode45',tol)).';
        end
        
        if nargout > 6 && calc_rp && calc_nr
             rho_gg_t(:,:,lp,j1,j2)   = tmp_gg_t + tmp_gg_cnj;
             rho_ee_t(:,:,lp,j1,j2) = tmp_ee_t+ tmp_ee_cnj;
        end 
        
        for j3 = 1:N
            if calc_nr

          R1tmp = V_eg_ee_left{j3}*tmp_ee_t;
          R4tmp = V_eg_gg_right{j3}*tmp_gg_t;
                    
            end
             if calc_rp
          
          R2tmp = V_eg_ee_left{j3}*tmp_ee_cnj;
          R3tmp = V_eg_gg_right{j3}*tmp_gg_cnj;
          
             end
            
        %    for j4 = j3   
        %in two colour experiment with om_1 = om_2 neq om_3 transitions 
        %with j3 neq j4 are heavily suppressed and are excluded           
             for j4 = 1:N   
              V_ge_sec = V_ge(:,:,j4);
            if calc_nr                  
                %first pure dipole factor
            ori_tmp1 = reshape(av1_g(j1,j2,j3,j4,:),[1,1,num_sets(1)]);
            %next the one which depends on |k_j|
            tmp2 = reshape(av2_g(j1,j2,j3,j4,:,:),size(av2_g,5),3);
            ori_tmp2 = reshape(e1(j1)*s_rp(1).*tmp2(:,1)+e1(j2)*s_rp(2).*tmp2(:,2)+...
                        e1(j3)*s_rp(3).*tmp2(:,3),[1,1,num_sets(2)]);
            ori_fct_rf = cat(3,ori_tmp1,ori_tmp2);
             %combine with the contributions, flip round t2 and t3    
            %this operation is a Liouville space inner product, i.e. a trace

            R_1 = mtimesx( mtimesx(V_ge_sec,R1tmp),ori_fct_rf); %t3_range by t2_range
            R_4 = mtimesx( mtimesx(V_ge_sec,R4tmp),ori_fct_rf);
            
            %also to preallocate R_6
            R_6 = zeros(size(R_1));
            end
            if calc_rp
            ori_tmp1 = reshape(av1_g(j1,j2,j3,j4,:),[1,1,num_sets(1)]);
            %next the one which depends on |k_j|
            tmp2 = reshape(av2_g(j1,j2,j3,j4,:,:),size(av2_g,5),3);
            ori_tmp2 = reshape(e1(j1)*s_nr(1).*tmp2(:,1)+e1(j2)*s_nr(2).*tmp2(:,2)+...
                        e1(j3)*s_nr(3).*tmp2(:,3),[1,1,num_sets(2)]);
%             ori_tmp2 = reshape((s_nr(1)+s_nr(2)).*tmp2(:,2).*om1_rng(lp)+...
%                         om3_rng*s_nr(3).*tmp2(:,3),[1,1,num_sets(2)]);
                                        
            ori_fct_nr = cat(3,ori_tmp1,ori_tmp2);                

            R_2 =  mtimesx(mtimesx(V_ge_sec,R2tmp),ori_fct_nr); %t3_range by t2_range
            R_3 =  mtimesx(mtimesx(V_ge_sec,R3tmp),ori_fct_nr);
            R_5 =  zeros(size(R_2));
            end       
                        
                for f3 = 1:N*(N-1)/2 %calculate excited state absorption cont
                    %for f4 = f3%1:N*(N-1)/2
                        %also excluded
                    for f4 = 1:N*(N-1)/2    
                    if calc_nr
                        R6tmp = -(V_fe_right{j3,f3}*tmp_ee_t); 
                 
             ori_tmp1 = reshape(av1_f(j1,j2,j3+(f3-1)*N,N*(f4-1)+j4,:),[1,1,num_sets(1)]);
            %next the one which depends on |k_j|
            tmp2 =reshape(av2_f(j1,j2,j3+(f3-1)*N,N*(f4-1)+j4,:,:),size(av2_f,5),3);
            ori_tmp2 = reshape(e1(j1)*s_rp(1).*tmp2(:,1)+e1(j2)*s_rp(2).*tmp2(:,2)+...
                        (e2(f3)-e1(j3))*s_rp(3).*tmp2(:,3),[1,1,num_sets(2)]);
            ori_fct_rp = cat(3,ori_tmp1,ori_tmp2); %combine together               
                   
            R6tmp2 = mtimesx( V_ef(:,:,j4,f4),R6tmp); %contract along dim 2,1
            %diff contribution for a range of dipole averages
            R_6 = R_6 + mtimesx( R6tmp2,ori_fct_rp); 
              %combine with the contributions 
                    end
                    if calc_rp
                        
                    R5tmp =  -(V_fe_right{j3,f3}*tmp_ee_cnj);  
                 
             ori_tmp1 = reshape(av1_f(j1,j2,j3+(f3-1)*N,N*(f4-1)+j4,:),[1,1,num_sets(1)]);
            %next the one which depends on |k_j|
            tmp2 =reshape(av2_f(j1,j2,j3+(f3-1)*N,N*(f4-1)+j4,:,:),size(av2_f,5),3);
            %assume |k_j| = c * omega_eg for transitions etc
            ori_tmp2 = reshape(e1(j1)*s_nr(1).*tmp2(:,1)+e1(j2)*s_nr(2).*tmp2(:,2)+...
                        (e2(f3)-e1(j3))*s_nr(3).*tmp2(:,3),[1,1,num_sets(2)]);
            ori_fct_nr = cat(3,ori_tmp1,ori_tmp2); %combine together
            
           R5tmp2 = mtimesx( V_ef(:,:,j4,f4),R5tmp);

            R_5 = R_5 + mtimesx( R5tmp2,ori_fct_nr);
                    end
                    end
                end
   %combine this frequency slices into the total response function             
               if calc_nr   
                   for lp2 = 1:size(R_2,3)
      R2(:,:,lp,lp2) = R2(:,:,lp,lp2) + R_2(:,:,lp2); 
      R3(:,:,lp,lp2) = R3(:,:,lp,lp2) + R_3(:,:,lp2); 
      R5(:,:,lp,lp2) = R5(:,:,lp,lp2) + R_5(:,:,lp2);
                   end
               end
               if calc_rp
                  for lp2 = 1:size(R_1,3)
      R1(:,:,lp,lp2) = R1(:,:,lp,lp2) + R_1(:,:,lp2); 
      R4(:,:,lp,lp2) = R4(:,:,lp,lp2) + R_4(:,:,lp2); 
      R6(:,:,lp,lp2) = R6(:,:,lp,lp2) + R_6(:,:,lp2);
                  end
               end
            end
        end
    end
    if t_verbose
    toc
    end
    end
end
end


%{
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
%}