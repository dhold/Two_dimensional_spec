function [R1,R2,R3,R4,R5,R6,rho_gg_t,rho_ee_t] = rp_cont_HEOM_NOA(sz_params,...
            rho_eg,rho_ge,V_set,V_ge_f,V_ef_f,prop_op_gg,prop_op_ee,...
            t2_range,calc_nr,calc_rp,H_trunc,tol,t_verbose) ;

% Calculated the re/nonrephasing contribution to a 2D signal given the 
% coherence part of the density matrix rho_ge(om), the Heisenberg picture 
% operators V(om_3) and the operators that propogate the system in the 
% ground and excited state.
% This version does no orientation averaging so the interaction operators
% are of the form \sum_j mu_xi_j.e_j V_xi_j, only the 2nd and 3rd ones
% need to be passed to this function in
%  V_set ={V_ge.e_2,V_ef.e_2,V_ge.e_3,V_ef.e_3)

%sz_params is a vector with the important length constants, note that code
%works (should) even if HL = 1 (not the HEOM, just a normal propogator)
N = sz_params(1); %number of sites
N2 = N*(N-1)/2; %number of double excited states\
HL = sz_params(2); %total number of aux density matricies (+1)
%lengths of the three different times
t1L = size(rho_eg,2);  t2L = length(t2_range);  t3L = size(V_ge_f,1);

%order of conjugation of pulses for each signal, 
%this is the sign in e^(+/- i(k dot r-omega t)), in order they interact
%s_rp = [-1,+1,+1]; %rephasing or photon echo
%s_nr = [+1,-1,+1]; %non rephasing signal

%pause on %nice to be able to stop the code in windows

if ~exist('tol','var'); tol =[1e-7,1e-5];
elseif isempty(tol);  tol =[1e-7,1e-5];
end %set tol if it isn't given

if ~exist('t_verbose','var')
    t_verbose = false; %don't output time for loops
end

if nargout > 6
   rho_gg_t = zeros(t2L,t1L,N,N);
   rho_ee_t  = zeros(N^2,t2L,t1L,N,N);
end
out_trunc_gg = H_trunc;
out_trunc_ee = N^2*H_trunc;

if ~calc_nr && ~calc_rp
   warning('nothing calculated')
   R1=[]; R2=[]; R3=[]; R4=[]; R5=[]; R6=[]; return
end
    
    if calc_nr  %preallocate nonrephasing
        R1 = complex(zeros(t3L,t2L,t1L));  R4 = R1; R6 = R4;
    else;  R1=[]; R4=[]; R6=[];   
    end
    if calc_rp   %preallocate rephasing
        R2 = complex(zeros(t3L,t2L,t1L));  R3 = R2; R5 = R3;
    else; R2 = []; R3=[]; R5 = [];  
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
iHL = speye(HL); %sparse identity matrix the size of the HEOM

   %this uses the arrow convention of "Advancing HEOM for efficient
   %evalution of ..." 1) is a right arrow (for op A) 2) is a left arrow
   %V_set ={V_ge.e_2,V_ef.e_2,V_eg.e_3,V_fe.e_3), V_set{2} only good for
   %coherence stuff

    V_ge = sparse(V_set{1});  V_eg = V_ge'; %second int
 
   V_eg_right = kron(iHL,kron(speye(N),V_eg));
   V_eg_left = kron(iHL,kron(V_eg.',speye(1))); 
   
   V_ge_right = kron(iHL,kron(speye(1),V_ge));
   V_ge_left = kron(iHL,kron((V_ge).',speye(N)));        

    V_eg = sparse(V_set{3}); % V_ge = V_eg'; %3rd int

   V_eg_gg_right = kron(iHL,kron(speye(1),V_eg));  
   V_eg_ee_left = kron(iHL,kron(V_eg.',speye(N)));   

   VV_fe = sparse(V_set{4});  %VV_ef = VV_fe';
   V_fe_right = kron(iHL, kron(speye(N),VV_fe));

    %de equations which deal with propogation
    rho_gg_prop = @(t,v) prop_op_gg*v;
    rho_ee_prop = @(t,v) prop_op_ee*v;

    if calc_nr  %preallocate nonrephasing
    tmp_gg = V_ge_right*rho_eg; 
    tmp_ee =  V_ge_left*rho_eg; 
    end
    if calc_rp   %preallocate rephasing
 % RP-> rho_ee = V_eg rho__ge , rho_gg = rho__ge V_eg 
    tmp_gg2 = V_eg_left*rho_ge;   %from conjugate path
    tmp_ee2 = V_eg_right*rho_ge; 
    end  

    for lp = 1:t1L %loop over t1L, WILL be painfully slow if t1L is large
%this loop could easily be parallelised with "parfor"
        if t_verbose;  tic; end
        
        if calc_nr
        tmp_gg_sec = tmp_gg(:,lp); %note this is a hole
        tmp_gg_t = (OD_wrapper(t2_range,rho_gg_prop,tmp_gg_sec,out_trunc_gg,'ode45',tol)).';      
        tmp_ee_sec = tmp_ee(:,lp); %excited state population
        tmp_ee_t = (OD_wrapper(t2_range,rho_ee_prop,tmp_ee_sec,out_trunc_ee,'ode45',tol)).';
        end

        if calc_rp
        tmp_gg_sec2 = tmp_gg2(:,lp);
        tmp_gg_cnj = (OD_wrapper(t2_range,rho_gg_prop,tmp_gg_sec2,out_trunc_gg,'ode45',tol)).';      

        tmp_ee_sec2 = tmp_ee2(:,lp);
        tmp_ee_cnj = (OD_wrapper(t2_range,rho_ee_prop,tmp_ee_sec2,out_trunc_ee,'ode45',tol)).';
        end
        
        if nargout > 6 && calc_rp && calc_nr
             rho_gg_t(:,lp)   = tmp_gg_t(1,:) + tmp_gg_cnj(1,:);
             rho_ee_t(:,:,lp) = tmp_ee_t(1:4,:) + tmp_ee_cnj(1:4,:);
        end 
        
            if calc_nr

          R1tmp = V_eg_ee_left*tmp_ee_t;
          R4tmp = V_eg_gg_right*tmp_gg_t;
          R6tmp = -V_fe_right*tmp_ee_t;    
                R1(:,:,lp) = V_ge_f*R1tmp; 
                R4(:,:,lp) = V_ge_f*R4tmp;
                R6(:,:,lp) = V_ef_f*R6tmp;

            end
             if calc_rp
          
          R2tmp = V_eg_ee_left*tmp_ee_cnj;
          R3tmp = V_eg_gg_right*tmp_gg_cnj;
          R5tmp = -V_fe_right*tmp_ee_cnj;  
                R2(:,:,lp) = V_ge_f*R2tmp; 
                R3(:,:,lp) = V_ge_f*R3tmp; 
                R5(:,:,lp) = V_ef_f*R5tmp;

             end
      
    if t_verbose;  toc
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