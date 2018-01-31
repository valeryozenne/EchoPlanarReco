
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
%
% Limitations:
%   only works with a single encoded space
%   fully sampled k-space (no partial fourier or undersampling)
%   one repetitions
%   doesn't handle averages, phases, segments and sets



clear all

addpath('/home/valery/Reseau/Valery/MatlabUnix/ismrm_sunrise_matlab-master/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading an existing file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = '/home/valery/Reseau/Imagerie/DICOM_DATA/2017-01-15_SMS/FID/00085_epi_RF_1dyn_SMS1_GRAPPA_FID10853.h5'
filename_noise = '/home/valery/Reseau/Imagerie/DICOM_DATA/2017-01-15_SMS/FID/00085_epi_RF_1dyn_SMS1_GRAPPA_NOISEFID10853.h5';
% filename_noise = '/home/valery/Reseau/Imagerie/DICOM_DATA/2017-01-15_SMS/FID/00087_epi_RF_1dyn_GRAPPA_NOISEFID10855.h5';
% filename = '/home/valery/Reseau/Imagerie/DICOM_DATA/2017-01-15_SMS/FID/00087_epi_RF_1dyn_GRAPPA_FID10855.h5';

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
% la technique consiste à annuler l’étape d’excitation RF dans le processus de formation du signal RMN. Ainsi, l’information électrique recueillie correspond exclusivement au bruit thermique
% la matrice de covariance de bruit se calcule à partir d’acquisitions de bruit.
% le calcul fait appel aux données brutes du bruit acquis par les différents éléments de réception. 
% En considérant que la matrice de données brutes est de taille, les éléments de la matrice de covariance s’exprimeront comme suit :
% La matrice de covariance comprend donc des termes qui traduisent les corrélations de bruit entre les éléments de réception (termes hors diagonale) ainsi que la variance du bruit capté par chaque récepteur (termes de la diagonale). L’on
% peut ainsi affirmer que les termes de la matrice de covariance de bruit montrent l’influence que les éléments de réception ont les uns par rapport aux autres" (extract from Anou SEWONU)




nbLines_noise=size(meas_noise.data,2);
% number of lines
str_msg=sprintf('le nombre TOTAL de lignes dans le fichier noise est %d \n', nbLines_noise); disp(str_msg);


%% calcul de la matrice de covariance
% pour commencer nous allons extraire les lignes correspondants au bruit
% en ouvrant l'ADC et acquiert le bruit electronique

% la première étape de traitement consiste à calculer la matrice de
% covariance du bruit
% cf fonction suivant du cours ismrm sunrise course 2013
% dmtx = ismrm_calculate_noise_decorrelation_mtx(noise(sp>0,:));

acq_noise_measurement =  find( (meas_noise.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT')) & ~(meas_noise.head.flagIsSet('ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA')) );
str_msg=sprintf('le nombre TOTAL de lignes dans noise  est %d \n', size(acq_noise_measurement,2)); disp(str_msg);

% nous allons utiliser 256 lignes

% allocation de la matrice de covariance du bruit
% la taille correspondant aux nombres d'antennes

dmtxTotal=zeros(nCoils,nCoils);

% pour chaque ligne, en tenant compte de l'oversampling, nous disposons
% d'une matrice à deux dimensions [enc_Nx, nCoils ]
% nous effections une permutation
% nous calculons le produit du bruit et de sa transposée
% nous sommons les données et comptons le nombre de points utilisés

Mtotal=0;

for p=1:size(acq_noise_measurement,2);
    
    noise_samples=meas_noise.data{acq_noise_measurement(p)};
    str_msg=sprintf('ligne %d  ,  matrix size [ %d , %d ] \n', p , size(noise_samples,1) , size(noise_samples,2));disp(str_msg);
    noise = reshape(noise_samples,numel(noise_samples)/size(noise_samples,length(size(noise_samples))), size(noise_samples,length(size(noise_samples))));
    noise = permute(noise,[2 1]);
    
    M = size(noise,2);
    % Performs generalized matrix multiplication.
    dmtx = (noise*noise');
    
    % nous comptons le nombre de points utilisés
    Mtotal=Mtotal+M;
    dmtxTotal=dmtxTotal+dmtx;
    
end

str_msg=sprintf('le nombre de points utilisés est %d   \n', Mtotal);disp(str_msg);

% normalisation
dmtxTotalNormalize = (1/(Mtotal-1))*dmtxTotal;


% figure matrix de covariance du bruit
figure(1)
imagesc(abs(dmtxTotal) , [0 1e-6]); title('matrix de covariance du bruit');


% Triangulation of matrices

% c'est important sinon cela ne correspond pas avec le gadgetron
% c'est pour annuler le terme sqrt(2) présent de le ismrm sunrise course

scale_factor=0.5;

dmtxTotalNormalize_temp = inv(chol(dmtxTotalNormalize,'lower'));
dmtxTotalNormalize_triangulation = dmtxTotalNormalize_temp*sqrt(2)*sqrt(scale_factor);





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
% clear D;


str_msg=sprintf('le nombre TOTAL de lignes dans le fichier image est %d \n', size(meas.data,2)); disp(str_msg);


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

noise_dwell_time_us_=meas_noise.head.sample_time_us(1);
acquisition_dwell_time_us_=meas.head.sample_time_us(1);
receiver_noise_bandwidth_=hdr.acquisitionSystemInformation.relativeReceiverNoiseBandwidth;
noise_bw_scale_factor_ = sqrt(2.0*acquisition_dwell_time_us_/noise_dwell_time_us_*receiver_noise_bandwidth_);

str_msg=sprintf('noise_dwell_time_us_ %f   \n', noise_dwell_time_us_);disp(str_msg);
str_msg=sprintf('acquisition_dwell_time_us_ %f   \n', acquisition_dwell_time_us_);disp(str_msg);
str_msg=sprintf('receiver_noise_bandwidth_ %f   \n', receiver_noise_bandwidth_);disp(str_msg);
str_msg=sprintf('noise_bw_scale_factor_ %f   \n', noise_bw_scale_factor_);disp(str_msg);


dmtxTotalNormalize_triangulation_bw=dmtxTotalNormalize_triangulation*noise_bw_scale_factor_;












measure_data_before_pre_whitening= zeros(enc_Nx*2, length(meas.data),  nCoils);

for p = 1:length(meas.data)
    
    measure_data_before_pre_whitening(:,p,:)=meas.data{p};
    
end


%%  apply the pre-whitening
in=measure_data_before_pre_whitening;
orig_size = size(in);
nelements = prod(orig_size)/nCoils;

disp(nelements)

in = permute(reshape(in,nelements,nCoils),[2 1]);
out = reshape(permute(dmtxTotalNormalize_triangulation_bw*in,[2 1]),orig_size);

data_pre_whitening=out;






% do the reco

contrast=1;
slice=1;
rep=1;
set=1;

% donne le nombre de ligne correspondant à ces parametres
acqs_image_only = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1))...
    & (meas.head.idx.set==(set-1))...
    & ~(meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) ...
    & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  );

str_msg=sprintf('le nombre de lignes ACQ_IS_IMAGE est  %d \n', size(acqs_image_only,2)); disp(str_msg);

% effectuons une reconstruction des données un peu brutale

kspace.raw = zeros(enc_Nx*2, enc_Ny,  nCoils);
kspace.pre_whitening = zeros(enc_Nx*2, enc_Ny,  nCoils);
mask.encoding= zeros(enc_Nx*2, enc_Ny);

for p = 1:length(acqs_image_only)
    ky = meas.head.idx.kspace_encode_step_1(acqs_image_only(p)) + 1;
    str_msg=sprintf('p %d  acqs_image_only(p)  %d  ky  %d \n', p, acqs_image_only(p), ky);  disp(str_msg);
    kspace.raw(:,ky,:) = meas.data{acqs_image_only(p)};
    kspace.pre_whitening(:,ky,:) = data_pre_whitening(:,acqs_image_only(p),:);
    mask.encoding(:,ky)=1;
end


figure(2)
subplot(2,2,1);imagesc(abs(ifft_2D(kspace.raw(:,:,1)))); colormap(gray);
subplot(2,2,2);imagesc(abs(ifft_2D(kspace.pre_whitening(:,:,1))));







%nous allons faire la correction EPI reconX

% il n'y a pas de fft a ce niveaux

contrast=1;
rep=1;

clear acq_reverse
acq_reverse=  find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.flagIsSet('ACQ_IS_REVERSE')) );

str_msg=sprintf('le nombre de lignes ACQ_IS_REVERSE est  %d \n', size(acq_reverse,2)); disp(str_msg);


clear acq_pas_reverse
acq_pas_reverse=  find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & ~(meas.head.flagIsSet('ACQ_IS_REVERSE')) );

str_msg=sprintf('le nombre de lignes ACQ_IS_PAS_REVERSE est  %d \n', size(acq_pas_reverse,2)); disp(str_msg);


reconx.encodeNx_  = hdr.encoding.encodedSpace.matrixSize.x;
reconx.encodeFOV_ = hdr.encoding.encodedSpace.fieldOfView_mm.x;
reconx.reconNx_   = hdr.encoding.reconSpace.matrixSize.x;
reconx.reconFOV_  = hdr.encoding.reconSpace.fieldOfView_mm.x;

temp_struct_long=hdr.encoding.trajectoryDescription.userParameterLong;
temp_struct_double=hdr.encoding.trajectoryDescription.userParameterDouble;

[ trajectoryPos_ , trajectoryNeg_, reconx ] = EPI_reconX_trapezoid_compute_trajectory( reconx, temp_struct_long , temp_struct_double);

[ Mpos_ ,   Mneg_  ]=EPI_reconX_trapezoid_compute_apply(   reconx , trajectoryPos_ , trajectoryNeg_, meas.head);


for p=1:1:length(acq_reverse)
    
    data.kspace1D.reconx(:,acq_reverse(p),:)= Mneg_*squeeze(data.kspace.pre_whitening(:,acq_reverse(p),:));
    
end

for p=1:1:length(acq_pas_reverse)
    
    data.kspace1D.reconx(:,acq_pas_reverse(p),:)= Mpos_*squeeze(data.kspace.pre_whitening(:,acq_pas_reverse(p),:));
    
end

data.kspace.reconx=fftshift(fft(fftshift(data.kspace1D.reconx,1),[],1),1);


kspace.reconx = zeros(enc_Nx, enc_Ny,  nCoils);

for p = 1:length(acqs_image_only)
    ky = meas.head.idx.kspace_encode_step_1(acqs_image_only(p)) + 1;
    str_msg=sprintf('p %d  acqs_image_only(p)  %d  ky  %d \n', p, acqs_image_only(p), ky);  disp(str_msg);    
    kspace.reconx(:,ky,:) = data.kspace.reconx(:,acqs_image_only(p),:);
   
end




figure(3)
subplot(2,2,1);imagesc(abs(ifft_2D(kspace.raw(:,:,1)))); colormap(gray);
subplot(2,2,2);imagesc(abs(ifft_2D(kspace.pre_whitening(:,:,1))));
subplot(2,2,3);imagesc(abs(ifft_2D(kspace.reconx(:,:,1))));






%% correction des ghosts de l'EPI

clear acq_phasecorr_data 
acq_phasecorr_data=  find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1)) ...
    & (meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) );
str_msg=sprintf('le nombre de lignes ACQ_IS_PHASECORR_DATA est  %d \n', size(acq_phasecorr_data,2)); disp(str_msg);

% il y a bien 2 * 3 lignes de correction pour l'EPI pour la coupe 1

clear acq_phasecorr_data_parallel_calibration
acq_phasecorr_data_parallel_calibration=  find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1)) ...
    & (meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA'))...
    & (meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION')));
str_msg=sprintf('le nombre de lignes ACQ_IS_PHASECORR_DATA_PARALLEL_CALIBRATION est  %d \n', size(acq_phasecorr_data_parallel_calibration,2)); disp(str_msg);

%
acq_parallel_calibration_only=  find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1)) ...
    & ~(meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA'))...
    & (meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION')));

str_msg=sprintf('le nombre de lignes ACQ_IS_PARALLEL_CALIBRATION_ONLY est  %d \n', size(acq_parallel_calibration_only,2)); disp(str_msg);

%
acq_phasecorr_data_image=  find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1)) ...
    & (meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA'))...
    & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION')));
str_msg=sprintf('le nombre de lignes ACQ_IS_PHASECORR_IMAGE est  %d \n', size(acq_phasecorr_data_image,2)); disp(str_msg);

%
acq_image_only=  find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1)) ...
    & ~(meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA'))...
    & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION')));

str_msg=sprintf('le nombre de lignes ACQ_IS_IMAGE_ONLY est  %d \n', size(acq_image_only,2)); disp(str_msg);


%% on va donc prendre les 3 premieres lignes

[ data_nav_phasecorr_calibration ] = read_and_store_navigator_data(meas,  acq_phasecorr_data_parallel_calibration , data_reconX);

figure()
c=1;
subplot(121); plot(abs(data_nav_phasecorr_calibration(:,1,c))); 
hold on ;  plot(abs(data_nav_phasecorr_calibration(:,2,c))); 
hold on ; plot(abs(data_nav_phasecorr_calibration(:,3,c)));
subplot(122); plot(angle(data_nav_phasecorr_calibration(:,1,c))); 
hold on ;  plot(angle(data_nav_phasecorr_calibration(:,2,c))); 
hold on ; plot(angle(data_nav_phasecorr_calibration(:,3,c)));



%% on va donc prendre juste la calibration

K1D_calibration = zeros(enc_Nx, enc_Ny,  nCoils);

for p = 1:length(acq_parallel_calibration_only)
    ky = meas.head.idx.kspace_encode_step_1(acq_parallel_calibration_only(p)) + 1;
    K1D_calibration(:,ky,:) = data_reconX(:,acq_parallel_calibration_only(p),:);
end

%% on va donc prendre juste l'image

K1D_image = zeros(enc_Nx, enc_Ny,  nCoils);

for p = 1:length(acq_image_only)
    ky = meas.head.idx.kspace_encode_step_1(acq_image_only(p)) + 1;
    K1D_image(:,ky,:) = data_reconX(:,acq_image_only(p),:);
end



%% calcul des coefficients pour la correction

readout=size(data_nav_phasecorr_calibration,1);
c=size(data_nav_phasecorr_calibration,3);

corrpos_=zeros( readout,1);
corrneg_=zeros( readout,1 );

[ slope, intercept, x ] = fit_slope( data_nav_phasecorr_calibration );

%% application des coefficients

tvec = slope*x + intercept;

corrpos_ = exp(complex(zeros( readout, 1), -0.5*tvec));
corrneg_ = exp(complex(zeros( readout, 1), +0.5*tvec));




clear acqs_phasecorr_calib_reverse clear acqs_phasecorr_calib_noreverse

acqs_phasecorr_calib_reverse = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1))...
    & (meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) ...
    & (meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  ...
    & (meas.head.flagIsSet('ACQ_IS_REVERSE'))   );


acqs_phasecorr_calib_noreverse = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1))...
    & (meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) ...
    & (meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  ...
    & ~(meas.head.flagIsSet('ACQ_IS_REVERSE'))   );

str_msg=sprintf('le nombre de lignes ACQ_IS_PARALLEL_CALIBRATION_REVERSE est  %d \n', size(acqs_phasecorr_calib_reverse,2)); disp(str_msg);
str_msg=sprintf('le nombre de lignes ACQ_IS_PARALLEL_CALIBRATION_NOREVERSE ACQ est  %d \n', size(acqs_phasecorr_calib_noreverse,2)); disp(str_msg);




data_nav_phasecorr_calibration_corrigee=zeros(size(data_nav_phasecorr_calibration));

% correction des lignes reverse
for p = 1:length(acqs_phasecorr_calib_reverse)    
    
      ky = meas.head.idx.kspace_encode_step_1(acqs_phasecorr_calib_reverse(p)) + 1
      str_msg=sprintf('ligne %d reverse   %d\n',p , acqs_phasecorr_calib_reverse(p)); disp(str_msg);
      
       for c=1:1:nCoils
            data_nav_phasecorr_calibration_corrigee(:,acqs_phasecorr_calib_reverse(p),c)=data_nav_phasecorr_calibration(:,acqs_phasecorr_calib_reverse(p),c).* corrneg_;
       end
end

% correction des lignes normale
for p = 1:length(acqs_phasecorr_calib_noreverse)    
    
    ky = meas.head.idx.kspace_encode_step_1(acqs_phasecorr_calib_noreverse(p)) + 1;
    str_msg=sprintf('ligne %d noreverse   %d\n',p , acqs_phasecorr_calib_noreverse(p)); disp(str_msg);
       for c=1:1:nCoils
            data_nav_phasecorr_calibration_corrigee(:,acqs_phasecorr_calib_noreverse(p),c)=data_nav_phasecorr_calibration(:,acqs_phasecorr_calib_noreverse(p),c).* corrpos_;
       end
end


figure(2)
c=1;
subplot(1,2,1)
plot(angle(data_nav_phasecorr_calibration(:,1,c)), 'r');
hold on;
plot(angle(data_nav_phasecorr_calibration(:,2,c)),'b');
hold on;
plot(angle(data_nav_phasecorr_calibration(:,3,c)),'g');

subplot(1,2,2)
plot(angle(data_nav_phasecorr_calibration_corrigee(:,1,c)), 'r');
hold on;
plot(angle(data_nav_phasecorr_calibration_corrigee(:,2,c)),'b');
hold on;
plot(angle(data_nav_phasecorr_calibration_corrigee(:,3,c)),'g');




%% on applique la correction à la calibration et à l'image

% attention ici nous n'avons pas pris la bonne correction, il faut prendre
% les 3 lignes de l'esapce de k précedent l'image et non la calibration acs

clear acqs_calib_reverse clear acqs_calib_noreverse

acqs_calib_reverse = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1))...
    & ~(meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) ...
    & (meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  ...
    & (meas.head.flagIsSet('ACQ_IS_REVERSE'))   );


acqs_calib_noreverse = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1))...
    & ~(meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) ...
    & (meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  ...
    & ~(meas.head.flagIsSet('ACQ_IS_REVERSE'))   );

str_msg=sprintf('le nombre de lignes ACQ_IS_PARALLEL_CALIBRATION_REVERSE est  %d \n', size(acqs_calib_reverse,2)); disp(str_msg);
str_msg=sprintf('le nombre de lignes ACQ_IS_PARALLEL_CALIBRATION_NOREVERSE ACQ est  %d \n', size(acqs_calib_noreverse,2)); disp(str_msg);



%% on va donc prendre juste la calibration

K1D_calibration_epi_corrected = zeros(enc_Nx, enc_Ny,  nCoils);

for p = 1:length(acqs_calib_reverse)    
    
    ky = meas.head.idx.kspace_encode_step_1(acqs_calib_reverse(p)) + 1
    
       for c=1:1:nCoils
            K1D_calibration_epi_corrected(:,ky,c)=K1D_calibration(:,ky,c).* corrneg_;
       end
end

for p = 1:length(acqs_calib_noreverse)    
    
    ky = meas.head.idx.kspace_encode_step_1(acqs_calib_noreverse(p)) + 1
    
       for c=1:1:nCoils
            K1D_calibration_epi_corrected(:,ky,c)=K1D_calibration(:,ky,c).* corrpos_;
       end
end



%% on va donc prendre les 3 premieres lignes de l'image

[ data_nav_phasecorr_image ] = read_and_store_navigator_data(meas,  acq_phasecorr_data_image , data_reconX);


%% calcul des coefficients pour la correction

readout=size(data_nav_phasecorr_image,1);
c=size(data_nav_phasecorr_image,3);

corrpos_=zeros( readout,1);
corrneg_=zeros( readout,1 );

[ slope, intercept, x ] = fit_slope( data_nav_phasecorr_image );

%% application des coefficients

tvec = slope*x + intercept;

corrpos_ = exp(complex(zeros( readout, 1), -0.5*tvec));
corrneg_ = exp(complex(zeros( readout, 1), +0.5*tvec));

%% on applique la correction à l'image

% attention ici nous n'avons pas pris la bonne correction, il faut prendre
% les 3 lignes de l'esapce de k précedent l'image et non la calibration acs

clear acqs_image_reverse clear acqs_image_noreverse

acqs_image_reverse = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1))...
    & ~(meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) ...
    & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  ...
    & (meas.head.flagIsSet('ACQ_IS_REVERSE'))   );


acqs_image_noreverse = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1))...
    & ~(meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) ...
    & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  ...
    & ~(meas.head.flagIsSet('ACQ_IS_REVERSE'))   );

str_msg=sprintf('le nombre de lignes ACQ_IS_IMAGE_REVERSE est  %d \n', size(acqs_image_reverse,2)); disp(str_msg);
str_msg=sprintf('le nombre de lignes ACQ_IS_IMAGE_NOREVERSE ACQ est  %d \n', size(acqs_image_noreverse,2)); disp(str_msg);



%% on va donc prendre juste l'image

K1D_image_epi_corrected = zeros(enc_Nx, enc_Ny,  nCoils);

for p = 1:length(acqs_image_reverse)    
    
    ky = meas.head.idx.kspace_encode_step_1(acqs_image_reverse(p)) + 1
    
       for c=1:1:nCoils
            K1D_image_epi_corrected(:,ky,c)=K1D_image(:,ky,c).* corrneg_;
       end
end

for p = 1:length(acqs_image_noreverse)    
    
    ky = meas.head.idx.kspace_encode_step_1(acqs_image_noreverse(p)) + 1
    
       for c=1:1:nCoils
            K1D_image_epi_corrected(:,ky,c)=K1D_image(:,ky,c).* corrpos_;
       end
end







% [ I_Image ] = ifft_2D(K_image);
% 
% [ I_Image_epi_corrected ] =  ifft_2D(K_image_epi_corrected);
% 
% 
% figure(10)
% for c = 1:nCoils
%     subplot(4,4,c); imagesc(abs(I_Image(:,:,c))); colormap(gray);
% end
% 
% figure(11)
% for c = 1:nCoils
%     subplot(4,4,c); imagesc(abs(I_Image_epi_corrected(:,:,c))); colormap(gray);
% end
% 
% figure(12)
% for c = 1:nCoils
%     subplot(4,4,c); imagesc(abs(I_Image_epi_corrected(:,:,c))-abs(I_Image(:,:,c))); colormap(gray);
% end

K_calibration_epi_corrected =   fftshift(fft(fftshift(K1D_calibration_epi_corrected,1),[],1),1);
K_calibration =   fftshift(fft(fftshift(K1D_calibration,1),[],1),1);


K_image_epi_corrected =   fftshift(fft(fftshift(K1D_image_epi_corrected,1),[],1),1);
K_image =   fftshift(fft(fftshift(K1D_image,1),[],1),1);


figure(10)
for c = 1:nCoils
    subplot(4,4,c); imagesc(abs(K_calibration_epi_corrected(:,:,c))); colormap(gray);
end

figure(11)
for c = 1:nCoils
    subplot(4,4,c); imagesc(abs(K_calibration(:,:,c))); colormap(gray);
end


K_for_grappa_epi_corrected=K_image_epi_corrected;
K_for_grappa=K_image;

for p = 1:length(acq_parallel_calibration_only)
    ky = meas.head.idx.kspace_encode_step_1(acq_parallel_calibration_only(p)) + 1;
    K_for_grappa_epi_corrected(:,ky,:)=K_calibration_epi_corrected(:,ky,:);
    K_for_grappa(:,ky,:)=K_calibration(:,ky,:);
end


figure(11)
for c = 1:nCoils
    subplot(4,4,c); imagesc(abs(K_for_grappa_epi_corrected(:,:,c))); colormap(gray);
end


sp = zeros(enc_Nx, enc_Ny);

for p = 1:length(acqs_image_only)
    ky = meas.head.idx.kspace_encode_step_1(acqs_image_only(p)) + 1;
    sp(:,ky) = 1;
    
end

for p = 1:length(acq_parallel_calibration_only)
    ky = meas.head.idx.kspace_encode_step_1(acq_parallel_calibration_only(p)) + 1;
    sp(:,ky) = sp(:,ky)+ 2;    
    sp(:,ky) = sp(:,ky+1)+ 2; 
end

figure()
imagesc(sp)



[img_grappa_coil_epi_corrected] = ismrm_cartesian_GRAPPA(K_for_grappa_epi_corrected,sp,2);
[img_grappa_coil] = ismrm_cartesian_GRAPPA(K_for_grappa,sp,2);

figure(11)
for c = 1:nCoils
    subplot(2,2,1); imagesc(abs(img_grappa_coil(:,:))); colormap(gray);
     subplot(2,2,2); imagesc(abs(img_grappa_coil_epi_corrected(:,:))); colormap(gray);
     
     subplot(2,2,3); imagesc(angle(img_grappa_coil(:,:))); colormap(gray);
     subplot(2,2,4); imagesc(angle(img_grappa_coil_epi_corrected(:,:))); colormap(gray);
end












%% partie PCA and coil compression (PCA coils gadgets)


%% pour effectuer la PCA rapidement sur le donnée , nous allons utliser
% uniquement les 100 premieres lignes de l'espace de k

% sur chacune de ces lignes nous allons garder les 16 point situés les plus
% au centre, par exemple, le point 64 etant le milieu du readout , nous
% gardons les points 57 à 72
% ce qui nous fait un total de 16*100 points * 30 antennes dans notre cas
% ensuite nous moyennons ces 1600 points pour chaque antenne et soustrayons
% la moyenne trouvé à chacun des points.


clear data_for_pca

for p = 1:length(acq_parallel_calibration_only)
    ky = meas.head.idx.kspace_encode_step_1(acq_parallel_calibration_only(p)) + 1;
    data_for_pca(:,p,:)=K_calibration_epi_corrected(:,ky,:);
    
end

% extraction des données dans un vecteur a part

size(data_for_pca)

sample_counter=0;
uncombined_channel=0;
samples_to_use=16;
data_offset=56;
channels=nCoils;
nbLinesForPca=length(acq_parallel_calibration_only);
samples_per_profile=128;
total_samples=nbLinesForPca*samples_to_use;
sample_counter=0;

clear A_ptr  tempo_v tempo_v2

means_ptr=zeros(nCoils,1);

for p = 1:length(acq_parallel_calibration_only)
    
    tempo_v=squeeze(data_for_pca(:,p,:));
    
    tempo_v2=reshape(tempo_v, [size(tempo_v,1)*size(tempo_v,2), 1] );
    
    size(tempo_v2)
    
    %ce nouveau vecteur a pour taille 128*30 =3840
    
    for  s =1:1:samples_to_use
        
        sample_counter=sample_counter+1;
        
        for  c = 1:1:nCoils
            
            indice_in=(c-1)*samples_per_profile + data_offset + s;
            
            A_ptr(sample_counter + (c-1)*total_samples,1)=tempo_v2(indice_in);
            
            means_ptr(c,1) = means_ptr(c,1) + tempo_v2(indice_in);
            
            %                 if (p==1 && c==11)%
            %                     str_msg=sprintf('%d c:  %1.5d %1.5d  \n',c ,   real(tempo_v2(indice_in)), imag(tempo_v2(indice_in)) ); disp(str_msg);%
            %                 end
                        str_msg=sprintf('Sample counter = %d/%d\n', sample_counter, total_samples); disp(str_msg);
        end
        
    end
end

size(A_ptr)

B_ptr=zeros(size(A_ptr));

% Subtraction de la moyenne
for  c = 1:1:nCoils
    
    for  s = 1:1:total_samples
        B_ptr(s + (c-1)*total_samples,1) =    A_ptr(s + (c-1)*total_samples,1) -means_ptr(c,1)/total_samples;
    end
    
end


%% Get eigenvectors

clear data_for_pca U S V

% en utilisant toutes les données
% data_for_pca = reshape(data_pre_whitening_remove,size(data_pre_whitening_remove,1)*size(data_pre_whitening_remove,2),nCoils);

% en utilisant uniquement les pixels au centre de l'espace des k
C_ptr=reshape(B_ptr, total_samples,nCoils);
data_for_pca = C_ptr;

[U S V]=svd(data_for_pca');
eigenvalues = S(:,1:nCoils).*S(:,1:nCoils);
eigenvectors = V(:,1:nCoils);

str_msg=sprintf('size eigenvalues   %d %d \n' , size(eigenvalues,1) , size(eigenvalues,2)); disp(str_msg);
str_msg=sprintf('size eigenvectors  %d %d \n' , size(eigenvectors,1) , size(eigenvectors,2)); disp(str_msg);
str_msg=sprintf('size U             %d %d \n' , size(U,1) , size(U,2)); disp(str_msg);

%the E is eigen value, the square of singular value

for c = 1:nCoils
    egienvalues_plot(c,1)=abs(eigenvalues(c,c));
end

figure()
plot(egienvalues_plot(:,1));


% si on utilise toute l'image on peut tracer les eigen vectors
%
% matrix_eigenvectors= reshape(eigenvectors,[size(data_pre_whitening_remove,1), size(data_pre_whitening_remove,2), nCoils]);
%
% image_matrix_eigenvectors=ifft_2D( matrix_eigenvectors );
%
% figure;
% for c = 1:nCoils
%     subplot(5,6,c); imagesc(abs(image_matrix_eigenvectors(:,:,c))); colormap(gray);
%     title(num2str(S(c,c)));
% end


% if (strcmp(do.check_agreement_with_gadgetron,'Y'))
%     
%     clear matrix_KT
%     
%     filename='/home/valery/Tempo/honDKLT_0_0.bin';
%     [ tempo_honDKLT ] = read_binary_complex( filename, nCoils*nCoils);
%     data_honDKLT=reshape(tempo_honDKLT, [nCoils , nCoils]);
%     matrix_KT=data_honDKLT;
%     
%     figure()
%     subplot(1,3,1); imagesc(abs(U));
%     subplot(1,3,2); imagesc(abs(matrix_KT));
%     subplot(1,3,3); imagesc(abs(U)-abs(matrix_KT));
%     
% end


K_calibration_epi_corrected_pca =    K_calibration_epi_corrected;
K_calibration_pca =  K_calibration;

K_image_epi_corrected_pca =    K_image_epi_corrected;
K_image_pca =    K_image;

for p = 1:size(K_calibration_epi_corrected,2)
    
    tempo_v=squeeze(K_calibration_epi_corrected(:,p,:));    
    new_tempo_v=tempo_v*U;    
    K_calibration_epi_corrected_pca(:,p,:)=new_tempo_v;
    
    
    tempo_v=squeeze(K_calibration(:,p,:));    
    new_tempo_v=tempo_v*U;    
    K_calibration_pca(:,p,:)=new_tempo_v;
    
end

for p = 1:size(K_calibration_epi_corrected,2)
    
    tempo_v=squeeze(K_image_epi_corrected(:,p,:));    
    new_tempo_v=tempo_v*U;    
    K_image_epi_corrected_pca(:,p,:)=new_tempo_v;
    
    
    tempo_v=squeeze(K_image(:,p,:));    
    new_tempo_v=tempo_v*U;    
    K_image_pca(:,p,:)=new_tempo_v;
    
end

I_image_epi_corrected_pca=ifft_2D(K_image_epi_corrected_pca);

I_image_pca=ifft_2D(K_image_pca);

figure(10)
for c = 1:nCoils
    subplot(4,4,c); imagesc(abs(I_image_epi_corrected_pca(:,:,c))); colormap(gray);
end


figure(11)
for c = 1:nCoils
    subplot(4,4,c); imagesc(abs(I_image_pca(:,:,c))); colormap(gray);
end


numberOfCoils=8;

K_for_grappa_epi_corrected_pca_reduced=K_image_epi_corrected_pca(:,:,1:numberOfCoils);
K_for_grappa_pca_reduced=K_image_pca(:,:,1:numberOfCoils);

K_calibration_epi_corrected_pca_reduced=K_calibration_epi_corrected_pca(:,:,1:numberOfCoils);
K_calibration_pca_pca_reduced=K_calibration_pca(:,:,1:numberOfCoils);


for p = 1:length(acq_parallel_calibration_only)
    ky = meas.head.idx.kspace_encode_step_1(acq_parallel_calibration_only(p)) + 1;
    K_for_grappa_epi_corrected_pca_reduced(:,ky,:)=K_calibration_epi_corrected_pca_reduced(:,ky,:);
    K_for_grappa_pca_reduced(:,ky,:)=K_calibration_pca_pca_reduced(:,ky,:);
end



[img_grappa_coil_epi_corrected_pca_reduced] = ismrm_cartesian_GRAPPA(K_for_grappa_epi_corrected_pca_reduced,sp,2);
[img_grappa_coil_pca_reduced] = ismrm_cartesian_GRAPPA(K_for_grappa_pca_reduced,sp,2);

figure(15)
for c = 1:nCoils
     subplot(2,2,1); imagesc(abs(img_grappa_coil_pca_reduced(:,:))); colormap(gray);
     subplot(2,2,2); imagesc(abs(img_grappa_coil_epi_corrected_pca_reduced(:,:))); colormap(gray);
     
     subplot(2,2,3); imagesc(angle(img_grappa_coil_pca_reduced(:,:))); colormap(gray);
     subplot(2,2,4); imagesc(angle(img_grappa_coil_epi_corrected_pca_reduced(:,:))); colormap(gray);
end






magnitudepre_whitening_remove_pca = zeros(enc_Nx,enc_Ny);
phase_pre_whitening_remove_pca = zeros(enc_Nx,enc_Ny);

sum_= zeros(enc_Nx,enc_Ny);
for c = 1:nCoils
    sum_ = sum_ + abs(im_pca_maison(:,:,c)).*exp(angle(im_pca_maison(:,:,c)).*1j);
end

magnitudepre_whitening_remove_pca = sqrt(sum(abs(im_pca_maison(:,:,:)).^2,3));
phase_pre_whitening_remove_pca = angle(sum_);



figure(10)
subplot(1,2,1);   imagesc(magnitudepre_whitening_remove_pca); colormap(gray)
subplot(1,2,2);   imagesc(phase_pre_whitening_remove_pca); colormap(gray)



figure(9)
for c = 1:nCoils
    subplot(5,6,c); imagesc(abs(im_pca_maison(:,:,c)), [0 1]); colormap(gray);
end

figure(10)
for c = 1:nCoils
    subplot(5,6,c); imagesc(angle(im_pca_maison(:,:,c)), [-pi pi]); colormap(gray);
end



%% check output pca


% %% Reconstruct images
% % Since the entire file is in memory we can use random access
% % Loop over repetitions, contrasts, slices
% reconImages = {};
% nimages = 0;
% for rep = 1:nReps
%     for contrast = 1:nContrasts
%         for slice = 1:nSlices
%             % Initialize the K-space storage array
%             K = zeros(2*enc_Nx, enc_Ny, enc_Nz, nCoils);
%             % Select the appropriate measurements from the data
%             acqs = find(  (meas.head.idx.contrast==(contrast-1)) ...
%                 & (meas.head.idx.repetition==(rep-1)) ...
%                 & (meas.head.idx.slice==(slice-1)));
%             for p = 1:length(acqs)
%                 ky = meas.head.idx.kspace_encode_step_1(acqs(p)) + 1;
%                 kz = meas.head.idx.kspace_encode_step_2(acqs(p)) + 1;
%                 K(:,ky,kz,:) = meas.data{acqs(p)};
%             end
%             % Reconstruct in x
%             K = fftshift(ifft(fftshift(K,1),[],1),1);
%             % Chop if needed
%             if (enc_Nx == rec_Nx)
%                 im = K;
%             else
%                 ind1 = floor((enc_Nx - rec_Nx)/2)+1;
%                 ind2 = floor((enc_Nx - rec_Nx)/2)+rec_Nx;
%                 im = K(ind1:ind2,:,:,:);
%             end
%             % Reconstruct in y then z
%             im = fftshift(ifft(fftshift(im,2),[],2),2);
%             if size(im,3)>1
%                 im = fftshift(ifft(fftshift(im,3),[],3),3);
%             end
%
%             % Combine SOS across coils
%             im = sqrt(sum(abs(im).^2,4));
%
%             % Append
%             nimages = nimages + 1;
%             reconImages{nimages} = im;
%         end
%     end
% end
%
% %% Display the first image
% figure
% colormap gray
% imagesc(reconImages{1}); axis image; axis off; colorbar;



