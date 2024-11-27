function [data_output  ] = apply_pca_basis(data_input, meas, acqs_image_all , U )

data_output=zeros(size(data_input));

for p=1:size(acqs_image_all,2)
    
    ky = meas.head.idx.kspace_encode_step_1(acqs_image_all(p)) + 1
    
    tempo_v=squeeze(data_input(:,ky,:));
    new_tempo_v=tempo_v*U;    
    data_output(:,ky,:)=new_tempo_v;
    
end

end

