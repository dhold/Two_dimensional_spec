%load(strcat('saved_data/plenio_test6_paramset_',num2str(2),'.mat'),'tmp_ee_t')


tmp_ee_t2=tmp_ee_t;
tmp_ee_t2(:,:,:,1,1) = tmp_ee_t2(:,:,:,1,1)*mu_ex(1,1)*mu_ex(1,1);
tmp_ee_t2(:,:,:,1,2) = tmp_ee_t2(:,:,:,1,2)*mu_ex(1,1)*mu_ex(2,1);
tmp_ee_t2(:,:,:,2,1) = tmp_ee_t2(:,:,:,2,1)*mu_ex(1,1)*mu_ex(2,1);
tmp_ee_t2(:,:,:,2,2) = tmp_ee_t2(:,:,:,2,2)*mu_ex(2,1)*mu_ex(2,1);
tmp_ee_t3 = tmp_ee_t2(:,:,:,1,1) + tmp_ee_t2(:,:,:,1,2)+tmp_ee_t2(:,:,:,2,1)+tmp_ee_t2(:,:,:,2,2);
tmp = tmp_ee_t3(:,:,25)+conj(tmp_ee_t3([1,3,2,4],:,25));
tmp = tmp/(tmp(1,1)+tmp(4,1));
tmp2 = reshape(tmp,2,2,length(t2_range));
figure
plot(t2_range,diagsum(mtimesx(tmp2,tmp2),1,2))