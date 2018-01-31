
%% Working with an existing ISMRMRD data set

% This is a simple example of how to reconstruct a single-shot EPI images from data
% acquired on a fully sampled cartesian grid
%
% Capabilities:
%   2D
%   use noise scans (pre-whitening)
%   remove oversampling in the readout direction
%   regridding of EPI traectory
%   ghost corretion using linear fit of phase diff between +/- readout
%   virtual coil using pca
%   coil reduction
%   magnitude reconstruction   %   
%
% Limitations:
%   only works with a single encoded space
%   no partial fourier
%   could handle repetitions averages, phases, segments and sets 
%   phase reconstruction is creepy

% Assuming
% 3 navigator lines at the beginning of acquisition
% a grappa 2 undersampling with calibration scan with 24 reference


% some arbitrary conventions 
% data is a structure that include data from the hdf5 without kspace reordering 
% the matrix dimensions are [R0 "shots" CHA], "shots" meaning the number of
% lines that has been acquired successively whatever the encoding  

% after noise prewhitening and regridding the data are reordered into a
% other structure called kspace, the matrix dimensions are [R0 E1 CHA]
% reconstruction is done for slice 1=1 , contrast =1, average =1 ,
% generalisation is easy using a loop

% so 
% data.kspace.raw is raw data
% data.kspace1D.raw is the raw data after ifft into the readout direction
% data.kspace1D.pre_whitening  is the data after  pre_whitening and with an ifft into the readout direction
% data.kspace1D.regridding  is the data after regridding and with an ifft into the readout direction

% kspace.image is the kspace of the image 
% kspace.calibration is the kspace of the calibration scan (24 lines usign Frappa 2 )   

% kspace.image_ghost_corrected is the kspace of the image after ghost correction 
% kspace.calibrationafter ghost correction  is the kspace of the
% calibration scan (24 lines using % Grappa 2 ) after ghost correction   




clear all

% download the ismrm_sunrise_matlab-master https://github.com/hansenms/ismrm_sunrise_matlab
% and change the link

addpath('/home/valery/Reseau/Valery/MatlabUnix/ismrm_sunrise_matlab-master/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading an existing file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% download the data and change the link below

% fully sample dataset
% available from the gadgetron in test integration using the python script
% get data.py 
% filename = '/home/valery/DICOM/FID/meas_MID517_nih_ep2d_bold_fa60_FID82077.h5';
% filename_noise = '/home/valery/DICOM/FID/meas_MID517_nih_ep2d_bold_fa60_NOISE82077.h5';
% figure_folder='figures_full_sampling/'; 


% undersampled dataset with grappa 2
% from my ftp server or send me an email
filename = '/home/valery/Reseau/Imagerie/DICOM_DATA/2017-01-15_SMS/FID/00085_epi_RF_1dyn_GRAPPA2_FID10853.h5'
filename_noise = '/home/valery/Reseau/Imagerie/DICOM_DATA/2017-01-15_SMS/FID/00085_epi_RF_1dyn_GRAPPA2_NOISEFID10853.h5';
figure_folder='figures_grappa_sampling/';  

mkdir(figure_folder) 

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
    nContrasts = hdr.encoding.encodingLimits.contrast.maximum + 1;
catch
    nContrasts = 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read all the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[meas , meas_noise ] = read_dset( dset_noise  , dset );

nbLines_totale=size(meas.data,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% apply noise pre whitening (NoiseAdjustGadget.cpp NoiseAdjustGadget.h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('apply noise pre whitening');

acq_noise_measurement =  find( (meas_noise.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT'))  ...
    & ~(meas_noise.head.flagIsSet('ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA')) );

str_msg=sprintf('le nombre TOTAL de lignes dans l acquisition de bruit  est %d \n', size(acq_noise_measurement,2)); disp(str_msg);
str_msg=sprintf('le nombre TOTAL de lignes dans le fichier image est %d \n', size(meas.data,2)); disp(str_msg);

oversampling=1;

if (size(acq_noise_measurement,2)>0)
    
    % si noise , extract_covariance_matrix
    [ triangular_covariance_matrix , dmtxTotalNormalize] = extract_covariance_matrix( meas_noise, acq_noise_measurement , nCoils );    
    
    [ data.kspace.pre_whitening ,  triangular_covariance_matrix_bw  ] = apply_pre_whitening( meas, meas_noise, hdr, triangular_covariance_matrix , oversampling );
    
else
    
    % else covariance matrice is the identity
    triangular_covariance_matrix=zeros(nCoils,nCoils);
    for c=1:nCoils
        triangular_covariance_matrix(c,c)=1;
    end    
    meas_noise.head.sample_time_us(1)=1;
    dmtxTotalNormalize=triangular_covariance_matrix;
    [ data.kspace.pre_whitening ,  triangular_covariance_matrix_bw  ] = apply_pre_whitening( meas, meas_noise, hdr, triangular_covariance_matrix , oversampling );    
    
end


figure(2)
subplot(1,2,2);imagesc(abs(triangular_covariance_matrix(:,:))); colormap(gray); title('triangular covariance matrix');
if (size(acq_noise_measurement,2)>0)
subplot(2,2,1);imagesc(abs(dmtxTotalNormalize(:,:)));   title('covariance matrix');
end
% subplot(2,2,3);imagesc(rsos(ifft_2D(kspace.reconx(:,:,:))));   title('after regridding');
print([figure_folder,'figure2'],'-dpng');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% apply gridding correction for ramp sampling (EPIReconXGadget.cpp EPIReconXGadget.h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('apply gridding correction');


% find the reverse line
clear acq_reverse
acq_reverse=  find((meas.head.flagIsSet('ACQ_IS_REVERSE')) );
str_msg=sprintf('le nombre de lignes ACQ_IS_REVERSE est  %d', size(acq_reverse,2)); disp(str_msg);

% find the non reverse line
clear acq_pas_reverse
acq_pas_reverse=  find( ~(meas.head.flagIsSet('ACQ_IS_REVERSE')) );
str_msg=sprintf('le nombre de lignes ACQ_IS_PAS_REVERSE est  %d', size(acq_pas_reverse,2)); disp(str_msg);

% reconx is a structure

reconx.encodeNx_  = hdr.encoding.encodedSpace.matrixSize.x;
reconx.encodeFOV_ = hdr.encoding.encodedSpace.fieldOfView_mm.x;
reconx.reconNx_   = hdr.encoding.reconSpace.matrixSize.x;
reconx.reconFOV_  = hdr.encoding.reconSpace.fieldOfView_mm.x;

temp_struct_long=hdr.encoding.trajectoryDescription.userParameterLong;
temp_struct_double=hdr.encoding.trajectoryDescription.userParameterDouble;

% calculate the trajectory
[ trajectoryPos_ , trajectoryNeg_, reconx ] = EPI_reconX_trapezoid_compute_trajectory( reconx, temp_struct_long , temp_struct_double);

disp('EPI_reconX_trapezoid_compute_trajectory');

[ Mpos_ ,   Mneg_  ] = EPI_reconX_trapezoid_compute_apply(   reconx , trajectoryPos_ , trajectoryNeg_, meas.head);

disp('EPI_reconX_trapezoid_compute_apply done');

data.kspace1D.reconx=zeros(size(data.kspace.pre_whitening,1)/2, size(data.kspace.pre_whitening,2), size(data.kspace.pre_whitening,3));

for p=1:1:length(acq_reverse)
    
    data.kspace1D.reconx(:,acq_reverse(p),:)= Mneg_*squeeze(data.kspace.pre_whitening(:,acq_reverse(p),:));
    
end

for p=1:1:length(acq_pas_reverse)
    
    data.kspace1D.reconx(:,acq_pas_reverse(p),:)= Mpos_*squeeze(data.kspace.pre_whitening(:,acq_pas_reverse(p),:));
    
end

data.kspace.reconx=fftshift(fft(fftshift(data.kspace1D.reconx,1),[],1),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% select and reorder the data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('select slice , contrat, average and reorder the data ');
% selection of the data

% you can create a loop here 
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

% copy and reorder the kspace according to the spatial encoding

kspace.raw = zeros(enc_Nx*2, enc_Ny,  nCoils);
kspace.pre_whitening = zeros(enc_Nx*2, enc_Ny,  nCoils);
kspace.reconx = zeros(enc_Nx, enc_Ny,  nCoils);

for p = 1:length(acqs_image_only)
    ky = meas.head.idx.kspace_encode_step_1(acqs_image_only(p)) + 1;
    
%     str_msg=sprintf('p %d  acqs_image_only(p)  %d  ky  %d \n', p, acqs_image_only(p), ky);  disp(str_msg);    
    
    kspace.raw(:,ky,:) = meas.data{acqs_image_only(p)};
    kspace.pre_whitening(:,ky,:) = data.kspace.pre_whitening(:,acqs_image_only(p),:);
    kspace.reconx(:,ky,:) = data.kspace.reconx(:,acqs_image_only(p),:);
   
end

% display the result

figure(3)
subplot(2,2,1);imagesc(rsos(ifft_2D(kspace.raw(:,:,:)))); colormap(gray); title('raw');
subplot(2,2,2);imagesc(rsos(ifft_2D(kspace.pre_whitening(:,:,:))));   title('after pre-whitening');
subplot(2,2,3);imagesc(rsos(ifft_2D(kspace.reconx(:,:,:))));   title('after regridding');
print([figure_folder,'figure3'],'-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% apply EPI ghost corrections (EPICorrGadget.cpp and EPICorrGadget.h with the linear fit option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('apply EPI ghost corrections');

% here, assuming a grappa 2 undersampling, we are selecting the lines/scans numbers corresponding to
% -1) EPI navigators for the calibration scan 
% -2) EPI navigators for the image
% -3) the calibration scan 
% -4) the image

% some steps could be removed 

%  just to check but useless
clear acq_phasecorr_data 
acq_phasecorr_data=  find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1)) ...
    & (meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) );
str_msg=sprintf('le nombre de lignes ACQ_IS_PHASECORR_DATA est  %d \n', size(acq_phasecorr_data,2)); disp(str_msg);

% il y a bien 2 * 3 lignes de correction pour l'EPI pour la coupe 1 (3 for
% the image + 3 for the calibration scans)


%-1)
clear acq_phasecorr_data_parallel_calibration
acq_phasecorr_data_parallel_calibration=  find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1)) ...
    & (meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA'))...
    & (meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION')));
str_msg=sprintf('le nombre de lignes ACQ_IS_PHASECORR_DATA_PARALLEL_CALIBRATION est  %d \n', size(acq_phasecorr_data_parallel_calibration,2)); disp(str_msg);

%-2)
acq_parallel_calibration_only=  find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1)) ...
    & ~(meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA'))...
    & (meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION')));
str_msg=sprintf('le nombre de lignes ACQ_IS_PARALLEL_CALIBRATION_ONLY est  %d \n', size(acq_parallel_calibration_only,2)); disp(str_msg);

%-3)
acq_phasecorr_data_image=  find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1)) ...
    & (meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA'))...
    & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION')));
str_msg=sprintf('le nombre de lignes ACQ_IS_PHASECORR_IMAGE est  %d \n', size(acq_phasecorr_data_image,2)); disp(str_msg);

%-4)
acq_image_only=  find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1)) ...
    & ~(meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA'))...
    & ~(meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION')));

str_msg=sprintf('le nombre de lignes ACQ_IS_IMAGE_ONLY est  %d \n', size(acq_image_only,2)); disp(str_msg);



if (hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1>1)


% on va donc prendre les 3 premieres lignes, read_and_store_navigator_data
% into a matrix [R0 number_of_navigator(3), CHA]

[ navigator.phasecorr_calibration ] = read_and_store_navigator_data(meas,  acq_phasecorr_data_parallel_calibration , data.kspace1D.reconx);


% display the navigator for channel 1, useless 

close(figure(4))
figure(4)
c=1;
subplot(121); plot(abs(navigator.phasecorr_calibration(:,1,c)));  title('navigator: magnitude');
hold on ;  plot(abs(navigator.phasecorr_calibration(:,2,c))); 
hold on ; plot(abs(navigator.phasecorr_calibration(:,3,c)));
subplot(122); plot(angle(navigator.phasecorr_calibration(:,1,c)));  title('navigator: phase');
hold on ;  plot(angle(navigator.phasecorr_calibration(:,2,c))); 
hold on ; plot(angle(navigator.phasecorr_calibration(:,3,c)));
print([figure_folder,'figure4'],'-dpng');


%% on va donc prendre juste la calibration (reordering the calibration scans into the kspace according to spatial encoding)

kspace1D.calibration = zeros(enc_Nx, enc_Ny,  nCoils);

for p = 1:length(acq_parallel_calibration_only)
    ky = meas.head.idx.kspace_encode_step_1(acq_parallel_calibration_only(p)) + 1;
    kspace1D.calibration(:,ky,:) = data.kspace1D.reconx(:,acq_parallel_calibration_only(p),:);
end

kspace.calibration=ifft_1D(kspace1D.calibration);

end



%% on va donc prendre juste l'image (reordering the image scans into the kspace according to spatial encoding)

kspace1D.image = zeros(enc_Nx, enc_Ny,  nCoils);

for p = 1:length(acq_image_only)
    ky = meas.head.idx.kspace_encode_step_1(acq_image_only(p)) + 1;
    kspace1D.image(:,ky,:) = data.kspace1D.reconx(:,acq_image_only(p),:);
end


kspace.image=ifft_1D(kspace1D.image);


%display the kspace and image
close(figure(5))
figure(5)
c=1;

subplot(222);imagesc(abs(kspace.image(:,:,c))); colormap(gray); title ('kspace image'); 
subplot(224);imagesc(rsos(ifft_2D(kspace.image(:,:,:))));  title ('image image')

if (hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1>1)
subplot(221);imagesc(abs(kspace.calibration(:,:,c)));  title ('kspace calib')
subplot(223);imagesc(rsos(ifft_2D(kspace.calibration(:,:,:))));   title ('image calibration')
end

print([figure_folder,'figure5'],'-dpng');
% here we want to estimate the phase shift to correct ghost niqist artefact 
% below the navigator scans of the calibration scans will be corrected, 
% in practice, it is useless, only the calibration scans must be corrected  


if (hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1>1)

%% calcul des coefficients pour la correction

readout=size(navigator.phasecorr_calibration,1);
c=size(navigator.phasecorr_calibration,3);

corrpos_=zeros( readout,1);
corrneg_=zeros( readout,1 );

% estimate slope from navigator scans
[ slope, intercept, x ] = fit_slope( navigator.phasecorr_calibration );

%% application des coefficients

tvec = slope*x + intercept;

corrpos_ = exp(complex(zeros( readout, 1), -0.5*tvec));
corrneg_ = exp(complex(zeros( readout, 1), +0.5*tvec));

% find the reverse and non reverse lines of phase corr scans in the
% calibration scans

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

% apply the ghost niquist correction after linear fit of phase difference
% between position and negative readout of the EPI

navigator.phasecorr_calibration_correted=zeros(size(navigator.phasecorr_calibration));

% correction des lignes reverse
for p = 1:length(acqs_phasecorr_calib_reverse)    
    
      ky = meas.head.idx.kspace_encode_step_1(acqs_phasecorr_calib_reverse(p)) + 1;      
%       str_msg=sprintf('ligne %d reverse   %d\n',p , acqs_phasecorr_calib_reverse(p)); disp(str_msg);
      
      for c=1:1:nCoils
            navigator.phasecorr_calibration_correted(:,acqs_phasecorr_calib_reverse(p),c)=navigator.phasecorr_calibration(:,acqs_phasecorr_calib_reverse(p),c).* corrneg_;
      end
end

% correction des lignes normales
for p = 1:length(acqs_phasecorr_calib_noreverse)    
    
    ky = meas.head.idx.kspace_encode_step_1(acqs_phasecorr_calib_noreverse(p)) + 1;
%     str_msg=sprintf('ligne %d noreverse   %d\n',p , acqs_phasecorr_calib_noreverse(p)); disp(str_msg);
    
       for c=1:1:nCoils
            navigator.phasecorr_calibration_correted(:,acqs_phasecorr_calib_noreverse(p),c)=navigator.phasecorr_calibration(:,acqs_phasecorr_calib_noreverse(p),c).* corrpos_;
       end
end

close(figure(6))
figure(6)
c=1;
subplot(1,2,1);  
plot(angle(navigator.phasecorr_calibration(:,1,c)), 'r'); title('no ghost niquist correction')
hold on;
plot(angle(navigator.phasecorr_calibration(:,2,c)),'b');
hold on;
plot(angle(navigator.phasecorr_calibration(:,3,c)),'g');

subplot(1,2,2)
plot(angle(navigator.phasecorr_calibration_correted(:,1,c)), 'r'); title('ghost niquist correction after linear fit of phase difference')
hold on;
plot(angle(navigator.phasecorr_calibration_correted(:,2,c)),'b');
hold on;
plot(angle(navigator.phasecorr_calibration_correted(:,3,c)),'g');
print([figure_folder,'figure6'],'-dpng');




% ok it is working now let's do the same for the calibration scans and next for the
% images
% on applique la correction à la calibration et à l'image

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

kspace1D.calibration_ghost_corrected = zeros(enc_Nx, enc_Ny,  nCoils);

for p = 1:length(acqs_calib_reverse)    
    
    ky = meas.head.idx.kspace_encode_step_1(acqs_calib_reverse(p)) + 1;
    
       for c=1:1:nCoils
            kspace1D.calibration_ghost_corrected(:,ky,c)=kspace1D.calibration(:,ky,c).* corrneg_;
       end
end

for p = 1:length(acqs_calib_noreverse)    
    
    ky = meas.head.idx.kspace_encode_step_1(acqs_calib_noreverse(p)) + 1;
    
       for c=1:1:nCoils
            kspace1D.calibration_ghost_corrected(:,ky,c)=kspace1D.calibration(:,ky,c).* corrpos_;
       end
end

kspace.calibration_ghost_corrected =  ifft_1D(kspace1D.calibration_ghost_corrected);

end


% attention ici nous n'avons pas pris la bonne correction, il faut prendre
% les 3 lignes de l'esapce de k précedent l'image et non la calibration acs
% warning , we need to redo the fit with the navigator scans from the image
% 

%% on va donc prendre les 3 premieres lignes de l'image

[ data_nav_phasecorr_image ] = read_and_store_navigator_data(meas,  acq_phasecorr_data_image , data.kspace1D.reconx);


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

kspace1D.image_ghost_corrected = zeros(enc_Nx, enc_Ny,  nCoils);

for p = 1:length(acqs_image_reverse)    
    
    ky = meas.head.idx.kspace_encode_step_1(acqs_image_reverse(p)) + 1;
    
       for c=1:1:nCoils
            kspace1D.image_ghost_corrected(:,ky,c)=kspace1D.image(:,ky,c).* corrneg_;
       end
end

for p = 1:length(acqs_image_noreverse)    
    
    ky = meas.head.idx.kspace_encode_step_1(acqs_image_noreverse(p)) + 1;
    
       for c=1:1:nCoils
            kspace1D.image_ghost_corrected(:,ky,c)=kspace1D.image(:,ky,c).* corrpos_;
       end
end




kspace.image_ghost_corrected =  ifft_1D(kspace1D.image_ghost_corrected);


close(figure(7))
figure(7)
c=1;

subplot(322);imagesc(abs(kspace.image(:,:,c)));  colormap(gray); title ('kspace image'); 
subplot(324);imagesc(rsos(ifft_2D(kspace.image(:,:,:))));  title ('image image')
subplot(326);imagesc(rsos(ifft_2D(kspace.image_ghost_corrected(:,:,:))));  title ('image image ghost corrected')

if (hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1>1)
    subplot(321);imagesc(abs(kspace.calibration(:,:,c))); colormap(gray); title ('kspace calib')
    subplot(323);imagesc(rsos(ifft_2D(kspace.calibration(:,:,:))));   title ('image calibration')
    subplot(325);imagesc(rsos(ifft_2D(kspace.calibration_ghost_corrected(:,:,:))));   title ('image calibration ghost corrected')
end
print([figure_folder,'figure7'],'-dpng');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% partie PCA and coil compression (PCACoilGadget.cpp PCACoilGadget.h  not really use anymore in the Generic Chain)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp('PCA and coil compression');

%% pour effectuer la PCA rapidement sur le donnée , nous allons utliser
% uniquement les 100 premieres lignes de l'espace de k

% sur chacune de ces lignes nous allons garder les 16 point situés les plus
% au centre, par exemple, le point 64 etant le milieu du readout , nous
% gardons les points 57 à 72
% ce qui nous fait un total de 16*100 points * 30 antennes dans notre cas
% ensuite nous moyennons ces 1600 points pour chaque antenne et soustrayons
% la moyenne trouvé à chacun des points.

if (hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1>1)
    
clear data_for_pca

for p = 1:length(acq_parallel_calibration_only)
    ky = meas.head.idx.kspace_encode_step_1(acq_parallel_calibration_only(p)) + 1;
    data_for_pca(:,p,:)=kspace.calibration_ghost_corrected(:,ky,:);
    
end

else
   
 clear data_for_pca

for p = round(enc_Ny/2)-12:1:round(enc_Ny/2)+11
    ky=p;
    data_for_pca(:,p-round(enc_Ny/2)+12+1,:)=kspace.image_ghost_corrected(:,ky,:);
    
end   
    
end


% extraction des données dans un vecteur a part

sample_counter=0;
uncombined_channel=0;
samples_to_use=16;
data_offset=enc_Nx/2-8;
channels=nCoils;
% nbLinesForPca=length(acq_parallel_calibration_only);
nbLinesForPca=size(data_for_pca,2);
samples_per_profile=enc_Ny;
total_samples=nbLinesForPca*samples_to_use;
sample_counter=0;

clear A_ptr  tempo_v tempo_v2

means_ptr=zeros(nCoils,1);

for p = 1:nbLinesForPca
        
    tempo_v=squeeze(data_for_pca(:,p,:));
    
    tempo_v2=reshape(tempo_v, [size(tempo_v,1)*size(tempo_v,2), 1] );        
    
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
%                         str_msg=sprintf('Sample counter = %d/%d\n', sample_counter, total_samples); disp(str_msg);
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

% display the eigenvector in ordre to find the cut off for the coil
% compression
close(figure(8))
figure(8)
plot(egienvalues_plot(:,1)); title('egienvalues after pca of the channels'); 
print([figure_folder,'figure8'],'-dpng');

% si on utilise toute l'image on peut tracer les eigen vectors

% matrix_eigenvectors= reshape(eigenvectors,[size(data_pre_whitening_remove,1), size(data_pre_whitening_remove,2), nCoils]);
% 
% image_matrix_eigenvectors=ifft_2D( matrix_eigenvectors );
% 
% figure;
% for c = 1:nCoils
%     subplot(5,6,c); imagesc(abs(image_matrix_eigenvectors(:,:,c))); colormap(gray);
%     title(num2str(S(c,c)));
% end



% apply the pca


if (hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1>1)
    
    kspace.calibration_ghost_corrected_pca =  kspace.calibration_ghost_corrected;
    kspace.calibration_pca =    kspace.calibration;
    
    for p = 1:size(kspace.image_ghost_corrected,2)
        
        tempo_v=squeeze(kspace.calibration_ghost_corrected(:,p,:));
        new_tempo_v=tempo_v*U;
        kspace.calibration_ghost_corrected_pca(:,p,:)=new_tempo_v;
        
        %useless
        tempo_v=squeeze(kspace.calibration(:,p,:));
        new_tempo_v=tempo_v*U;
        kspace.calibration_pca(:,p,:)=new_tempo_v;
        
    end
    
end


kspace.image_ghost_corrected_pca =    kspace.image_ghost_corrected;
kspace.image_pca =    kspace.image;


for p = 1:size(kspace.image_ghost_corrected,2)
    
    tempo_v=squeeze(kspace.image_ghost_corrected(:,p,:));    
    new_tempo_v=tempo_v*U;    
    kspace.image_ghost_corrected_pca(:,p,:)=new_tempo_v;
    
    %useless 
    tempo_v=squeeze(kspace.image(:,p,:));    
    new_tempo_v=tempo_v*U;    
    kspace.image_pca(:,p,:)=new_tempo_v;
    
end


% display the result

figure(10)
for c = 1:nCoils
    subplot(4,4,c); imagesc(abs(ifft_2D(kspace.image_ghost_corrected(:,:,c)))); colormap(gray); title(['real channels image ', num2str(c)]);
end
print([figure_folder,'figure10'],'-dpng');

figure(11)
for c = 1:nCoils
    subplot(4,4,c); imagesc(abs(ifft_2D(kspace.image_ghost_corrected_pca(:,:,c)))); colormap(gray); title(['virtual channels image ', num2str(c)]);
end
print([figure_folder,'figure11'],'-dpng');

if (hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1>1)

figure(12)
for c = 1:nCoils
    subplot(4,4,c); imagesc(abs(ifft_2D(kspace.calibration_ghost_corrected(:,:,c)))); colormap(gray); title(['real channels calibration ', num2str(c)]);
end
print([figure_folder,'figure12'],'-dpng');

figure(13)
for c = 1:nCoils
    subplot(4,4,c); imagesc(abs(ifft_2D(kspace.calibration_ghost_corrected_pca(:,:,c)))); colormap(gray); title(['virtual channels calibration ', num2str(c)]);
end
print([figure_folder,'figure13'],'-dpng');
end

%%%%%%%%%%%%%%%%%%%%%%%%
%% coil compression (CoilReductionGadget.cpp CoilReductionGadget.h)
%%%%%%%%%%%%%%%%%%%%%%%%%%

number_of_coil_to_keep=4; 
% take the egienvalues data after pca and find the virtual channel for
% which its eigenvector is 1% of the eigenvector of the first virtual channel 


kspace.image_ghost_corrected_reduced=kspace.image_ghost_corrected_pca(:,:,1:number_of_coil_to_keep);
kspace.image_reduced=kspace.image_pca(:,:,1:number_of_coil_to_keep);

if (hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1>1)

kspace.calibration_ghost_corrected_reduced=kspace.calibration_ghost_corrected_pca(:,:,1:number_of_coil_to_keep);
kspace.calibration_reduced=kspace.calibration_pca(:,:,1:number_of_coil_to_keep);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simple 2D ifft recontruction
%%%%%%%%%%%%%%%%%%%%%%%%%%

if (hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1==1)
    
figure(17)
subplot(2,2,1); imagesc(rsos(ifft_2D(kspace.image_ghost_corrected_reduced))); colormap(gray); title('all raw data')
subplot(2,2,2); imagesc(rsos(ifft_2D(kspace.image_ghost_corrected))); colormap(gray); title(['using only ',num2str(number_of_coil_to_keep),' virtual channels over ', num2str(channels)]);
subplot(2,2,3); imagesc(rsos(ifft_2D(kspace.image_ghost_corrected_reduced))-rsos(ifft_2D(kspace.image_ghost_corrected))); colormap(gray); title('diff');
print([figure_folder,'figure17'],'-dpng');
end


if (hdr.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1>1)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% grappa recontruction
%%%%%%%%%%%%%%%%%%%%%%%%%%

% WARNING this the grappa for the sunrise course, there is some differences with the Grappa of the gadgetron 
% WARNING calibration shoudl be used independently to estimate the kernel
% and then apply to the folded image
% nevertheless, here, there is no motion or major changes in whatever 
% so it "works"
 

% defined the mask for the grappa

sp = zeros(enc_Nx, enc_Ny);
sp_image = zeros(enc_Nx, enc_Ny);
sp_calib = zeros(enc_Nx, enc_Ny);

for p = 1:length(acqs_image_only)
    ky = meas.head.idx.kspace_encode_step_1(acqs_image_only(p)) + 1;
    sp(:,ky) = 1;
    sp_image(:,ky) = 1;
    
end

for p = 1:length(acq_parallel_calibration_only)
    ky = meas.head.idx.kspace_encode_step_1(acq_parallel_calibration_only(p)) + 1;
    sp(:,ky) = sp(:,ky)+ 2;    
    sp(:,ky) = sp(:,ky+1)+ 2; 
    
    sp_calib(:,ky) = sp_calib(:,ky)+ 2;
    sp_calib(:,ky) = sp_calib(:,ky+1)+ 2;
    
end

% display the encoding pattern
close(figure(15))
figure(15)
imagesc(sp)
print([figure_folder,'figure15'],'-dpng');
% merge the calibration scan and the image scan

acqs_parallel = find(  (meas.head.idx.contrast==(contrast-1)) ...
    & (meas.head.idx.repetition==(rep-1)) ...
    & (meas.head.idx.slice==(slice-1))...
    & (meas.head.idx.set==(set-1))...
    & ~(meas.head.flagIsSet('ACQ_IS_PHASECORR_DATA')) ...
    & (meas.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION'))  );

% put the missing calibration scan into the kspace
kspace_for_grappa=kspace.image_ghost_corrected_reduced;
for p = 1:length(acqs_parallel)
    ky = meas.head.idx.kspace_encode_step_1(acqs_parallel(p)) + 1;
    kspace_for_grappa(:,ky,:)=kspace.calibration_ghost_corrected_reduced(:,ky,:);

end

%display the encoding pattern and the corresponding kspace 
close(figure(16))
figure(16)
c=1;
subplot(231);imagesc(abs(kspace.calibration_ghost_corrected_reduced(:,:,c)), [0 15]); colormap(gray); title ('kspace calib');
subplot(232);imagesc(abs(kspace.image_ghost_corrected_reduced(:,:,c)), [0 15]);  title ('kspace image');
subplot(233);imagesc(abs(kspace_for_grappa(:,:,c)), [0 15]);   title ('kspace image + calib');
subplot(234);imagesc(sp_image); colormap(gray); title ('mask kspace calib');
subplot(235);imagesc(sp_calib);  title ('mask kspace image');
subplot(236);imagesc(sp);   title ('mask  kspace image + calib');
print([figure_folder,'figure16'],'-dpng');
% do the grappa calculation
[image.ghost_corrected_reduced] = ismrm_cartesian_GRAPPA(kspace_for_grappa,sp,2);

clear kspace_for_grappa
% repeat without coil compression
kspace_for_grappa=kspace.image_ghost_corrected;
for p = 1:length(acqs_parallel)
    ky = meas.head.idx.kspace_encode_step_1(acqs_parallel(p)) + 1;
    kspace_for_grappa(:,ky,:)=kspace.calibration_ghost_corrected(:,ky,:);

end

% do the grappa calculation
[image.ghost_corrected] = ismrm_cartesian_GRAPPA(kspace_for_grappa,sp,2);

figure(17)
subplot(2,2,1); imagesc(abs(image.ghost_corrected)); colormap(gray); title('all raw data')
subplot(2,2,2); imagesc(abs(image.ghost_corrected_reduced)); colormap(gray); title(['using only ',num2str(number_of_coil_to_keep),' virtual channels over ', num2str(channels)]);
subplot(2,2,3); imagesc(abs(image.ghost_corrected)-abs(image.ghost_corrected_reduced)); colormap(gray); title('diff');
  print([figure_folder,'figure17'],'-dpng');  
% display the creepy phase images
% figure(16)
% subplot(2,2,1); imagesc(angle(image.ghost_corrected)); colormap(gray); title('all raw data')
% subplot(2,2,2); imagesc(angle(image.ghost_corrected_reduced)); colormap(gray); title(['using only ',num2str(number_of_coil_to_keep),' virtual channels over ', num2str(size(kspace_for_grappa,3))]);
% subplot(2,2,3); imagesc(angle(image.ghost_corrected)-abs(image.ghost_corrected_reduced)); colormap(gray); title('diff');
     


end
