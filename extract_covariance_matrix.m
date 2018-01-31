function [ dmtxTotalNormalize_triangulation , dmtxTotalNormalize, dmtxTotal] = extract_covariance_matrix( meas_noise , acq_noise_measurement , nCoils , debug)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here







% nous allons utiliser 256 lignes

% allocation de la matrice de covariance du bruit, la taille correspondant aux nombres d'antennes



% pour chaque ligne, en tenant compte de l'oversampling, nous disposons
% d'une matrice à deux dimensions [enc_Nx, nCoils ]
% nous effections une permutation
% nous calculons le produit du bruit et de sa transposée
% nou sommons les données et comptons le nombre de points utilisés

dmtxTotal=zeros(nCoils,nCoils);

Mtotal=0;

for p=1:size(acq_noise_measurement,2);
    
    noise_samples=double(meas_noise.data{acq_noise_measurement(p)});
%     str_msg=sprintf('ligne %d  ,  matrix size [ %d , %d ] \n', p , size(noise_samples,1) , size(noise_samples,2));disp(str_msg);
    noise = reshape(noise_samples,numel(noise_samples)/size(noise_samples,length(size(noise_samples))), size(noise_samples,length(size(noise_samples))));
%     noise = permute(noise,[2 1]);
    
    M = size(noise,1);
    
    % Performs generalized matrix multiplication.
    dmtx = (noise'*noise);
       
    
    % nous comptons le nombre de points utilisés pour la normalisation
    % future
    Mtotal=Mtotal+M;
    dmtxTotal=dmtxTotal+dmtx;
    
end

str_msg=sprintf('le nombre de points utilisés est %d  ', Mtotal);disp(str_msg);

% normalisation
dmtxTotalNormalize = (1/(Mtotal-1))*dmtxTotal;

% Triangulation of matrices

% le facteur est important sinon cela ne correspond pas avec le Gadgetron
% TODO trouver une meilleur raison

scale_factor=0.5;

dmtxTotalNormalize_temp = tril(inv(chol(dmtxTotalNormalize,'lower')));
dmtxTotalNormalize_triangulation = dmtxTotalNormalize_temp*sqrt(2)*sqrt(scale_factor);

% figure matrix de covariance du bruit
% figure(1)
% subplot(121); imagesc(abs(dmtxTotalNormalize)); title('matrix de covariance du bruit');
% subplot(122); imagesc(abs(dmtxTotalNormalize_triangulation)); title('matrix de covariance du bruit');



end

