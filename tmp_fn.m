LL = sz1*Nv;
R_red_op = zeros(LL^2);

for a = 1:sz1
    for b = 1:sz1

    temp = sparse(squeeze(R_red(a,b,:,:)));
        for av = 1:Nv
            for bv = 1:Nv
                temp2 = sparse(av,bv,1,Nv,Nv); %this vibrational element
                temp3 = kron(temp,temp2); %combine vib and elec
                
    lp = ((b-1)*Nv + bv-1)*sz1*Nv + (a-1)*Nv + av; %this element is mapped to

    R_red_op(lp,:) = reshape(temp3,[LL^2,1]);
            end
        end
    end
end