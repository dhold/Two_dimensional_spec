function [Supop,Lindblad_op] = redfield_lindblad_calc2(H_e,H_f,H_vib,...
            fock_space_rep,M_prj,Beta,Drude_modes,BO_modes)

        N = length(H_e);
if ~isnan(H_f)   %pass nan to take only 1st ex manifold
sz1 =  1+length(H_e) + length(H_f);
to_pass = blkdiag(1,H_f);
else
 sz1 =  length(H_e); to_pass = [];
end
sz2 =  length(H_vib);

R_red = redfield_calc_2(H_e,to_pass,...
            fock_space_rep,Beta,Drude_modes,[],[],[],[],true);
%Supop = R_red; Lindblad_op = []; return        
[Lindblad_op] = Lindblad_op_gen(Beta,BO_modes,[],[],sz1,M_prj);

LL = sqrt(length(Lindblad_op));

R_red_op = zeros(LL^2); 

temp = zeros(LL);
cnt1 = 1; cnt2 = 1;
cnt1vib =0; cnt2vib =1; 
%warning('this bit seems wrong')
%cnt_test = 0;
for lp = 1:LL^2
    
    %iterate along one vibrational lvl
    cnt1vib = cnt1vib+1;
    if cnt1vib > sz2  %if beyond this level reset to one and add another to cnt1
        cnt1vib=1; cnt1 = cnt1+1;
        if cnt1 > sz1
            cnt1 = 1;  cnt2vib = cnt2vib+1;
            if cnt2vib > sz2
                cnt2vib = 1; cnt2 = cnt2+1;
            end
        end
    end
    temp = temp*0;
   for a = 1:sz1 %loop over electronic degrees of freedom
        for b = 1:sz1 
            temp((a-1)*sz2+cnt1vib,(b-1)*sz2+cnt2vib) = R_red(cnt1,cnt2,a,b);
        end
   end
%cnt_test = cnt_test + sum(sum(temp~=0));
    R_red_op(lp,:) = reshape(temp,1,sz1^2*sz2^2 ).';

end
%cnt_test
%numel(R_red_op)

R_red_op = sparse(R_red_op);
if nargout == 1
Supop =  -R_red_op + sparse(Lindblad_op);
end

sz_full = sz1*sz2;
if ~isnan(H_f) %cut up into sections
    tmpg = zeros(sz_full); tmpg(1:sz2,1:sz2)=1;
	tmpg = logical(reshape(tmpg,sz_full^2,1)); %picks this out

%reduced operators acting only on ground excited coherences, these don't
%mix to p_gg or p_ee'

    tmpeg = zeros(sz_full); tmpeg(1:sz2,sz2+1:sz2*(N+1))=1;
    %tmpeg(sz2+1:sz2*(N+1),1:sz2)=1;
    tmpeg = logical(reshape(tmpeg,sz_full^2,1));
    
%reduced operators acting only on 1st ex state manifold, 

    tmpe  = zeros(sz_full); 
    tmpe(sz2+1:sz2*(N+1),sz2+1:sz2*(N+1)) = 1;
	tmpe = logical(reshape(tmpe,sz_full^2,1));
    
%reduced operator acting on gs-double ex coherence

%     tmpgf  = zeros(sz_full); tmpgf(1:sz2,sz2*(N+1)+1:end)=1;
%     %tmpgf(sz2*(N+1)+1:end,1:sz2)=1;
%     tmpgf = logical(reshape( tmpgf,sz_full^2,1));
    
%reduced operators acting only on 1st-2nd ex state manifold, 

    tmpef  = zeros(sz_full); tmpef(sz2+1:sz2*(N+1),sz2*(N+1)+1:end)=1; %upper diag
   % tmpef(sz2*(N+1)+1:end,sz2+1:sz2*(N+1))=1; %lower diag
	tmpef = logical(reshape(tmpef,sz_full^2,1));

    if nargout == 1
Supop = {Supop(tmpg,tmpg),Supop(tmpeg,tmpeg),Supop(tmpe,tmpe),Supop(tmpef,tmpef)};
    else
Supop = {R_red_op(tmpg,tmpg),R_red_op(tmpeg,tmpeg),R_red_op(tmpe,tmpe),R_red_op(tmpef,tmpef)};
Lindblad_op = {Lindblad_op(tmpg,tmpg),Lindblad_op(tmpeg,tmpeg),Lindblad_op(tmpe,tmpe),Lindblad_op(tmpef,tmpef)};         
    end
end