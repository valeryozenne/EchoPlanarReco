
%% Working with an existing ISMRMRD data set

% This is a simple example of how to reconstruct images from data
% acquired on a fully sampled cartesian grid
%
% Capabilities:
%   2D
%   use noise scans (pre-whitening)
%   remove oversampling in the readout direction
%   virtual coil using pca
%   coil reduction
%   magnitude and phase reconstruction
%   read data output from gadgetron
%
% Limitations:
%   only works with a single encoded space
%   fully sampled k-space (no partial fourier or undersampling)
%   one repetitions
%   doesn't handle averages, phases, segments and sets
%
%

% We first create a data set using the example program like this:
%   ismrmrd_generate_cartesian_shepp_logan -r 5 -C -o shepp-logan.h5
% This will produce the file shepp-logan.h5 containing an ISMRMRD
% dataset sampled evenly on the k-space grid -128:0.5:127.5 x -128:127
% (i.e. oversampled by a factor of 2 in the readout direction)
% with 8 coils, 5 repetitions and a noise level of 0.5
% with a noise calibration scan at the beginning


clear all


addpath('/home/valeryozenne/mount/Valery/MatlabBAART/Toolboxes/Generic_functions/');

[ str_user ] = get_PC_name();
[ str_network_imagerie, str_network_perso, str_network_dicom ] = get_network_name( str_user );

%addpath('../Utility');
%addpath('../Noise');

addpath(['/home/',str_user,str_network_perso ,'/MatlabUnix/ismrm_sunrise_matlab-master/']);
addpath('/usr/local/share/ismrmrd/matlab/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading an existing file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filename = '/home/valeryozenne/mount/Imagerie/For_Valery/TP_reco/FID/meas_MID00026_FID04585_gre_with_phase.h5';
filename_noise = '/home/valeryozenne/mount/Imagerie/For_Valery/TP_reco/FID/meas_MID00026_NOISE04585_gre_with_phase.h5';

if exist(filename_noise, 'file')
    dset_noise = ismrmrd.Dataset(filename_noise, 'dataset');
else
    error(['File ' filename_noise ' does not exist.  Please generate it.'])
end


if exist(filename, 'file')
    dset = ismrmrd.Dataset(filename, 'dataset');
else
    error(['File ' filename ' does not exist.  Please generate it.'])
end

check_agreement_with_gadgetron='N';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read some fields from the XML header %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We need to check if optional fields exists before trying to read them

hdr_noise = ismrmrd.xml.deserialize(dset_noise.readxml);
hdr = ismrmrd.xml.deserialize(dset.readxml);


%% Encoding and reconstruction information
% Matrix size
enc_Nx = hdr.encoding.encodedSpace.matrixSize.x;
enc_Ny = hdr.encoding.encodedSpace.matrixSize.y;
enc_Nz = hdr.encoding.encodedSpace.matrixSize.z;
rec_Nx = hdr.encoding.reconSpace.matrixSize.x;
rec_Ny = hdr.encoding.reconSpace.matrixSize.y;
rec_Nz = hdr.encoding.reconSpace.matrixSize.z;

% Field of View
enc_FOVx = hdr.encoding.encodedSpace.fieldOfView_mm.x;
enc_FOVy = hdr.encoding.encodedSpace.fieldOfView_mm.y;
enc_FOVz = hdr.encoding.encodedSpace.fieldOfView_mm.z;
rec_FOVx = hdr.encoding.reconSpace.fieldOfView_mm.x;
rec_FOVy = hdr.encoding.reconSpace.fieldOfView_mm.y;
rec_FOVz = hdr.encoding.reconSpace.fieldOfView_mm.z;

% Number of slices, coils, repetitions, contrasts etc.
% We have to wrap the following in a try/catch because a valid xml header may
% not have an entry for some of the parameters

try
    nSlices = hdr.encoding.encodingLimits.slice.maximum + 1;
catch
    nSlices = 1;
end

try
    nCoils = hdr.acquisitionSystemInformation.receiverChannels;
catch
    nCoils = 1;
end

try
    nReps = hdr.encoding.encodingLimits.repetition.maximum + 1;
catch
    nReps = 1;
end

try
    nContrasts = hdr.encoding.encodingLimits.contrast.maximum + 1 + 1;
catch
    nContrasts = 1;
end


% TODO add the other possibilites

%% Read all the data
% Reading can be done one acquisition (or chunk) at a time,
% but this is much faster for data sets that fit into RAM.
D_noise = dset_noise.readAcquisition();
D = dset.readAcquisition();





%% Take noise scans for pre-whitening example

% Find the first non-noise scan
% This is how to check if a flag is set in the acquisition header
isNoise = D_noise.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');

% toutes les données sont du bruit

firstScan = find(isNoise==1,1,'first');
if firstScan > 1
    noise = D_noise.select(1:firstScan-1);
else
    noise = [];
end

meas_noise  = D_noise.select(firstScan:D_noise.getNumber);
clear D_noise;



%% quelques vérifications

if (hdr.encoding.encodedSpace.matrixSize.x~=hdr.encoding.reconSpace.matrixSize.x)
    
    disp(' presence d oversampling \n');
    
end


%% pre-whitening (NoiseAdjustGadget.cpp NoiseAdjustGadget.h)

nbLines_noise=size(meas_noise.data,2);

%nb de lignes
str_msg=sprintf('le nombre TOTAL de lignes dans le fichier noise est %d \n', nbLines_noise); disp(str_msg);


%% calcul de la matrice de covariance
% pour commencer nous allons extraire les lignes correspondants au bruit
% ouvre l'ADC et acquiert le bruit electronique

% la première étape de traitement consiste à calculer la matrice de
% covariance du bruit
% cf fonction suivant du cours ismrm sunrise course 2013
% dmtx = ismrm_calculate_noise_decorrelation_mtx(noise(sp>0,:));

acq_noise_measurement =  find( (meas_noise.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT')) & ~(meas_noise.head.flagIsSet('ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA')) );
str_msg=sprintf('le nombre TOTAL de lignes dans noise  est %d \n', size(acq_noise_measurement,2)); disp(str_msg);

% nous allons utiliser 256 lignes

% allocation de la matrice de covariance du bruit
% la taille correspondant aux nombres d'antennes



% pour chaque ligne, en tenant compte de l'oversampling, nous disposons
% d'une matrice à deux dimensions [enc_Nx, nCoils ]
% nous effections une permutation
% nous calculons le produit du bruit et de sa transposée
% nous sommons les données et comptons le nombre de points utilisés

[ mat.covariance_matrix.triangular ,  mat.covariance_matrix.raw_normalise, mat.covariance_matrix.raw] = extract_covariance_matrix( meas_noise, acq_noise_measurement , nCoils );


%% Ignore noise scans
% TODO add a pre-whitening example
% Find the first non-noise scan
% This is how to check if a flag is set in the acquisition header
isNoise = D.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
firstScan = find(isNoise==0,1,'first');
if firstScan > 1
    noise = D.select(1:firstScan-1);
else
    noise = [];
end

meas  = D.select(firstScan:D.getNumber);
clear D;

size(meas.data)

% If the goal of the noise decorrelation is to create data with noise SD = 1.0 in
%  real and imaginary channels, the scale_factor should be:
%
%     scale_factor = (T_acq_dwell/T_noise_dwell)*NoiseReceiverBandwidthRatio
%
%  Where T_acq_dwell is the sampling dwell time in the acquisition data being
%  decorrelated and T_noise_dwell is the dwell time in the noise sample
%  acquistion (they may not be the same). The NoiseReceiverBandwithRatio is
%  a number (typically smaller than 1) that takes into account that the
%  noise acquisition may not have flat frequency response.
%
%  see Kellman P, McVeigh ER. Magn Reson Med. 2005 Dec;54(6):1439-47. Erratum in: Magn Reson Med. 2007 Jul;58(1):211-2.

% les valeurs suivantes sont necessaires

% noise_dwell_time_us_=meas_noise.head.sample_time_us(1);
% acquisition_dwell_time_us_=meas.head.sample_time_us(1);
% receiver_noise_bandwidth_=hdr.acquisitionSystemInformation.relativeReceiverNoiseBandwidth;
% noise_bw_scale_factor_ = sqrt(2.0*acquisition_dwell_time_us_/noise_dwell_time_us_*receiver_noise_bandwidth_);
%
% str_msg=sprintf('noise_dwell_time_us_ %f   \n', noise_dwell_time_us_);disp(str_msg);
% str_msg=sprintf('acquisition_dwell_time_us_ %f   \n', acquisition_dwell_time_us_);disp(str_msg);
% str_msg=sprintf('receiver_noise_bandwidth_ %f   \n', receiver_noise_bandwidth_);disp(str_msg);
% str_msg=sprintf('noise_bw_scale_factor_ %f   \n', noise_bw_scale_factor_);disp(str_msg);



[ mat.data.pre_whitening ,  mat.covariance_matrix.triangular_bw  , mat.data.before_pre_whitening] = apply_pre_whitening( meas, meas_noise, hdr, mat.covariance_matrix.triangular, 0 );



%% do the reco


contrast=1;
slice=1;
rep=1;
set=1;


% donne le nombre de ligne correspondant à ces parametres
acqs_image_all = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1))...
    & (meas.head.idx.set==(set-1))   );

str_msg=sprintf('le nombre de lignes ACQ_IS_ ALL est  %d \n', size(acqs_image_all,2)); disp(str_msg);



% on passe les donnees en kspace

mat.kspace.raw = zeros(enc_Nx, enc_Ny,  nCoils);
mat.kspace.pre_whitening = zeros(enc_Nx, enc_Ny,  nCoils);
mat.sampling = zeros(enc_Nx, enc_Ny);


for p = 1:length(acqs_image_all)
    
    ky = meas.head.idx.kspace_encode_step_1(acqs_image_all(p)) + 1;
    str_msg=sprintf('p %d  acqs_image_all(p)  %d  ky  %d', p, acqs_image_all(p), ky);  disp(str_msg);
    
    mat.kspace.raw(:,ky,:) = meas.data{acqs_image_all(p)};
    mat.kspace.pre_whitening(:,ky,:) = mat.data.pre_whitening(:,acqs_image_all(p),:);
    mat.sampling(:,ky)=1;
    
end


c=1;
figure(21)

subplot (2,3,1) ; imagesc(abs(mat.kspace.raw(:,:,c)));colormap(gray);
subplot (2,3,2) ; imagesc(angle(mat.kspace.raw(:,:,c)), [-pi pi]);
subplot (2,3,4) ; imagesc(abs(mat.kspace.pre_whitening(:,:,c)));
subplot (2,3,5) ; imagesc(angle(mat.kspace.pre_whitening(:,:,c)), [-pi pi]);
    
    
%% reco after remove oversampling pour montrer ce que cela donne

% gt.image.pre_whitening= ifft_2D(gt.kspace.pre_whitening);
mat.image.pre_whitening= ifft_2D(mat.kspace.pre_whitening*sqrt(2));

mat.image.remove_oversampling=zeros(enc_Nx/2,enc_Ny, nCoils);

mat.image.remove_oversampling=  mat.image.pre_whitening((enc_Nx/4)+1:enc_Nx*3/4,:,:);

mat.kspace.remove_oversampling=fft_2D(mat.image.remove_oversampling);



c=1;
figure(26)
subplot (2,3,2) ; imagesc(abs(mat.kspace.remove_oversampling(:,:,c)));
subplot (2,3,5) ; imagesc(angle(mat.kspace.remove_oversampling(:,:,c)), [-pi pi]);
subplot (2,3,1) ; imagesc(abs(mat.image.remove_oversampling(:,:,c))); colormap(gray);
subplot (2,3,4) ; imagesc(angle(mat.image.remove_oversampling(:,:,c) ), [-pi pi]);


%% partie PCA and coil compression (PCA coils gadgets)


% pour effectuer la PCA rapidement sur le donnée , nous allons utliser
% uniquement les 100 premieres lignes de l'espace de k

% sur chacune de ces lignes nous allons garder les 16 point situés les plus
% au centre, par exemple, le point 64 etant le milieu du readout , nous
% gardons les points 57 à 72
% ce qui nous fait un total de 16*100 points * 30 antennes dans notre cas
% ensuite nous moyennons ces 1600 points pour chaque antenne et soustrayons
% la moyenne trouvé à chacun des points.

% todo : attention la fonction calculate_pca_basis differe un peu dans la
% version generic grappa

% calcul de la matrice de transfert

[ U , S, V , eigenvalues, eigenvalues_plot, nbLinesMaximumForPCA] = calculate_pca_basis_generic(mat.kspace.remove_oversampling, meas,  acqs_image_all );

% on applique la nouvelle base

[ mat.kspace.pca ] = apply_pca_basis(mat.kspace.remove_oversampling, meas, acqs_image_all , U );

% reconstruction des données sans grappa après pca juste pour montrer ce
% que cela donne

[ mat.image.pca ] = ifft_2D( mat.kspace.pca );




figure(31)
for c=1:nCoils
subplot(6,6,c); imagesc(abs(ifft_2D(mat.kspace.remove_oversampling(:,:,c)))); colormap(gray);
end


figure(32)
for c=1:nCoils
subplot(6,6,c); imagesc(abs(ifft_2D(mat.kspace.pca(:,:,c)))); colormap(gray);
end




mat.image.magnitude.sos = zeros(enc_Nx/2,enc_Ny);
mat.image.phase.sos = zeros(enc_Nx/2,enc_Ny);

sum_= zeros(enc_Nx/2,enc_Ny);
for c = 1:nCoils
    sum_ = sum_ + abs(mat.image.pca (:,:,c)).*exp(angle(mat.image.pca (:,:,c)).*1j);
end

mat.image.magnitude.sos = sqrt(sum(abs(mat.image.pca (:,:,:)).^2,3));
mat.image.phase.sos = angle(sum_);

figure(10)
subplot(1,2,1);   imagesc(mat.image.magnitude.sos); colormap(gray)
subplot(1,2,2);   imagesc(mat.image.phase.sos); colormap(gray)




figure(9)
for c = 1:nCoils
    subplot(5,6,c); imagesc(abs(mat.image.pca(:,:,c)), [0 1]); colormap(gray);
end

figure(10)
for c = 1:nCoils
    subplot(5,6,c); imagesc(angle(mat.image.pca(:,:,c)), [-pi pi]); colormap(gray);
end





%% nous allons maintenant attaquer les gadgets GrappaGadget GrappaUnimixingGadget

% le kspace est strictement identique entre maltab et gadgetron



sp = zeros(enc_Nx/2, enc_Ny);
sp2 = zeros(enc_Nx/2, enc_Ny);

for p = 1:length(acqs_image_only)
    ky = meas.head.idx.kspace_encode_step_1(acqs_image_only(p)) + 1;
    sp(:,ky) = 1;
    sp2(:,ky) = 1;
end

for p = 1:length(acqs_parallel_calibration)
    ky = meas.head.idx.kspace_encode_step_1(acqs_parallel_calibration(p)) + 1;
    sp(:,ky) = sp(:,ky)+ 2;
    sp2(:,ky) = sp2(:,ky)+ 2;
    sp2(:,ky+1) = sp2(:,ky+1)+ 2;
end




nb_ok_coil_to_keep=8;

[mat.image.grappa.wash_pd] = ismrm_cartesian_GRAPPA(mat.kspace.pca(:,:,1:nb_ok_coil_to_keep),sp,2);

[mat.image.grappa.wash] = ismrm_cartesian_GRAPPA(mat.kspace.pca(:,:,1:nb_ok_coil_to_keep),sp,2);

[mat.image.grappa.wash_kellman] = ismrm_cartesian_GRAPPA_val(mat.kspace.pca(:,:,1:nb_ok_coil_to_keep),sp,2, [], [] , 1);

[mat.image.grappa.wash_python] = ismrm_cartesian_GRAPPA_val(mat.kspace.pca(:,:,1:nb_ok_coil_to_keep),sp,2, [], [] , 2);

[mat.image.grappa.inati_iter_python] = ismrm_cartesian_GRAPPA_val(mat.kspace.pca(:,:,1:nb_ok_coil_to_keep),sp,2, [], [] , 3);

[mat.image.grappa.inati_gt] = ismrm_cartesian_GRAPPA_val(mat.kspace.pca(:,:,1:nb_ok_coil_to_keep),sp,2, [], [] , 4);

[mat.image.grappa.inati_matlab] = ismrm_cartesian_GRAPPA_val(mat.kspace.pca(:,:,1:nb_ok_coil_to_keep),sp,2, [], [] , 3);



figure()
subplot(2,4,1);  imagesc(abs(mat.image.grappa.wash_kellman));
subplot(2,4,2);  imagesc(angle(mat.image.grappa.wash_kellman)); colormap(gray);
subplot(2,4,3);  imagesc(abs(mat.image.grappa.wash_python));
subplot(2,4,4);  imagesc(angle(mat.image.grappa.wash_python)); colormap(gray);
subplot(2,4,5);  imagesc(abs(mat.image.grappa.inati_matlab));
subplot(2,4,6);  imagesc(angle(mat.image.grappa.inati_matlab)); colormap(gray);


filename_magnitude_out5='/home/valery/Dev/Data/cpu_grappa_flash_m.real';
filename_phase_out5='/home/valery/Dev/Data/cpu_grappa_flash_p.real';

mat.image.magnitude.thruth=read_binary_array(filename_magnitude_out5, 2, [128 128]);
mat.image.phase.thruth=read_binary_array(filename_phase_out5, 2, [128 128]);

max(mat.image.magnitude.thruth(:))
min(mat.image.magnitude.thruth(:))

max(abs(mat.image.grappa.inati_matlab(:)))
min(abs(mat.image.grappa.inati_matlab(:)))


figure(20)
subplot(2,3,1); imagesc(mat.image.magnitude.thruth); colormap(gray); title('truth');
subplot(2,3,2); imagesc(abs(mat.image.grappa.inati_matlab)); colormap(gray); title('test');
subplot(2,3,3); imagesc(mat.image.magnitude.thruth-abs(mat.image.grappa.inati_matlab)); colormap(gray);
subplot(2,3,4); imagesc(image_phase_truth); colormap(gray); title('truth');
subplot(2,3,5); imagesc(angle(mat.image.grappa.inati_matlab)); colormap(gray); title('test');
subplot(2,3,6); imagesc(image_phase_truth-angle(mat.image.grappa.inati_matlab)); colormap(gray);










%% check data avant coil estimation inati
if (strcmp(check_agreement_with_gadgetron,'Y'))
    
%     for c = 1:16
%         
%         str_c=sprintf('%d',c-1);
%         filename=['/home/valery/Tempo/target_acs_',str_c,'.bin'];
%         [ tempo_gt ] = read_binary_complex( filename, enc_Nx/2*enc_Ny);
%         target_acs(:,:,c)=reshape(tempo_gt, [enc_Nx/2 , enc_Ny]);
%         
%     end
    
%     for c = 1:16
%         
%         str_c=sprintf('%d',c-1);
%         filename=['/home/valery/Tempo/complex_im_',str_c,'.bin'];
%         [ tempo_gt ] = read_binary_complex( filename, enc_Nx/2*enc_Ny);
%         complex_im(:,:,c)=reshape(tempo_gt, [enc_Nx/2 , enc_Ny]);
%         
%     end
    
%     for c = 1:16
%         
%         str_c=sprintf('%d',c-1);
%         filename=['/home/valery/Tempo/coil_map_',str_c,'.bin'];
%         [ tempo_gt ] = read_binary_complex( filename, enc_Nx/2*enc_Ny);
%         coil_map(:,:,c)=reshape(tempo_gt, [enc_Nx/2 , enc_Ny]);
%         
%     end
    
    
%     figure(10);
%     for c = 1:1
%         subplot(4,3,1);  imagesc(abs(target_acs(:,:,c)));  colormap(gray)
%         subplot(4,3,2);  imagesc(abs(mat.kspace.pca_klt(:,:,c)));
%         subplot(4,3,3);  imagesc(abs(target_acs(:,:,c))-abs(mat.kspace.pca_klt(:,:,c)));
%         subplot(4,3,4);  imagesc(abs(ifft_2D(target_acs(:,:,c))));
%         subplot(4,3,5);  imagesc(abs(ifft_2D(mat.kspace.pca_klt(:,:,c))));
%         subplot(4,3,6);  imagesc(abs(ifft_2D(target_acs(:,:,c)))-abs(ifft_2D(mat.kspace.pca_klt(:,:,c))));
%         subplot(4,3,7);  imagesc(angle(ifft_2D(target_acs(:,:,c))));
%         subplot(4,3,8);  imagesc(angle(ifft_2D(mat.kspace.pca_klt(:,:,c))));
%         subplot(4,3,9);  imagesc(angle(ifft_2D(target_acs(:,:,c)))-angle(ifft_2D(mat.kspace.pca_klt(:,:,c))));
%         subplot(4,3,10);  imagesc(angle(complex_im(:,:,c)));
%         subplot(4,3,11);  imagesc(angle(ifft_2D(target_acs(:,:,c))));
%         subplot(4,3,12);  imagesc(angle(complex_im(:,:,c))-angle(ifft_2D(target_acs(:,:,c))));
%         pause(0.2)
%     end
        
end



%  [ coilMap_matlab ] = coil_map_study_2d_Inati( complex_im, 5, 3 );


% filename=['/home/valery/Tempo/D_','0_0','.bin'];
% [ tempo_gt ] = read_binary_complex( filename, 25*16);
% D_gt(:,:,c)=reshape(tempo_gt, [25 , 16]);


% 
% figure(10);
% for c = 1:1
%     subplot(2,3,1);  imagesc(abs(coil_map(:,:,c)));  colormap(gray)
%     subplot(2,3,2);  imagesc(abs(coilMap_matlab(:,:,c)));
%     subplot(2,3,3);  imagesc(abs(coil_map(:,:,c))-abs(coilMap_matlab(:,:,c)));
%     subplot(2,3,4);  imagesc(angle(coil_map(:,:,c)));
%     subplot(2,3,5);  imagesc(angle(coilMap_matlab(:,:,c)));
%     subplot(2,3,6);  imagesc(angle(coil_map(:,:,c))-angle(coilMap_matlab(:,:,c)));
% end








