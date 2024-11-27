function [ U , S, V , eigenvalues, eigenvalues_plot, nbLinesMaximumForPCA] = calculate_pca_basis_generic(data_input, meas, acqs_image )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% nbLinesMaximumForPCA=100;

nbLinesMaximumForPCA=2; %length(acqs_image);


% uncombined_channel=0;
samples_to_use=16;
data_offset=56;
nCoils=size(data_input,3);
samples_per_profile=size(data_input,1);
total_samples=nbLinesMaximumForPCA*samples_to_use;
sample_counter=0;

clear A_ptr  tempo_v tempo_v2

means_ptr=zeros(nCoils,1);

compteur=1;

data_for_pca=zeros(size(data_input,1), nbLinesMaximumForPCA , size(data_input,3));

for p=1:nbLinesMaximumForPCA
    
    ky = meas.head.idx.kspace_encode_step_1(acqs_image(p)) + 1;
    
    data_for_pca(:,compteur,:)=data_input(:,ky,:);
    compteur=compteur+1;
end


for p=1:nbLinesMaximumForPCA
    
    tempo_v=squeeze(data_for_pca(:,p,:));  % 128 *1* 20  on prend la première ligne
    
    tempo_v2=reshape(tempo_v, [size(tempo_v,1)*size(tempo_v,2), 1] );  %un vecteur 2560 
    
    %ce nouveau vecteur a pour taille 128*30 =3840
    
    for  s =1:1:samples_to_use
        
        sample_counter=sample_counter+1;
        
        for  c = 1:1:nCoils
            
            indice_in=(c-1)*samples_per_profile + data_offset + s;
            
            A_ptr(sample_counter + (c-1)*total_samples,1)=tempo_v2(indice_in);
            
            means_ptr(c,1) = means_ptr(c,1) + tempo_v2(indice_in);            
          
        end
        
    end
end

clear B_ptr

% Subtraction de la moyenne
for  c = 1:1:nCoils
    
    for  s = 1:1:total_samples
        B_ptr(s + (c-1)*total_samples,1) =    A_ptr(s + (c-1)*total_samples,1) -means_ptr(c,1)/total_samples;
    end
    
end




%% Get eigenvectors

clear data_for_pca U S V eigenvalues eigenvectors eigenvalues_plot eigenv_plot

% en utilisant toutes les données
% data_for_pca = reshape(data_pre_whitening_remove,size(data_pre_whitening_remove,1)*size(data_pre_whitening_remove,2),nCoils);

% en utilisant uniquement les pixels au centre de l'espace des k
C_ptr=reshape(B_ptr, total_samples,nCoils);
size(C_ptr)
data_for_pca = C_ptr;
size(data_for_pca')

[U, S, V]=svd(data_for_pca');
eigenvalues = S(:,1:nCoils).*S(:,1:nCoils);
eigenv = S(:,1:nCoils);
eigenvectors = V(:,1:nCoils);

str_msg=sprintf('size eigenvalues   %d %d' , size(eigenvalues,1) , size(eigenvalues,2)); disp(str_msg);
str_msg=sprintf('size eigenvectors  %d %d' , size(eigenvectors,1) , size(eigenvectors,2)); disp(str_msg);
str_msg=sprintf('size U             %d %d' , size(U,1) , size(U,2)); disp(str_msg);

%the E is eigen value, the square of singular value

for c = 1:nCoils
    eigenvalues_plot(c,1)=eigenvalues(c,c);
    eigenv_plot(c,1)=eigenv(c,c);
end


return
