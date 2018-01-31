function [ out  , triangular_covariance_matrix_bw , data_before_pre_whitening ] = apply_pre_whitening( meas, meas_noise, hdr, triangular_covariance_matrix , oversampling)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


% les valeurs suivantes sont necessaires

enc_Nx = hdr.encoding.encodedSpace.matrixSize.x;
nCoils = hdr.acquisitionSystemInformation.receiverChannels;

noise_dwell_time_us_=meas_noise.head.sample_time_us(1);
acquisition_dwell_time_us_=meas.head.sample_time_us(1);
receiver_noise_bandwidth_=hdr.acquisitionSystemInformation.relativeReceiverNoiseBandwidth;
noise_bw_scale_factor_ = sqrt(2.0*acquisition_dwell_time_us_/noise_dwell_time_us_*receiver_noise_bandwidth_);

str_msg=sprintf('enc_Nx %f   ', enc_Nx);disp(str_msg);
str_msg=sprintf('nCoils %f   ', nCoils);disp(str_msg);
str_msg=sprintf('noise_dwell_time_us_ %f   ', noise_dwell_time_us_);disp(str_msg);
str_msg=sprintf('acquisition_dwell_time_us_ %f   ', acquisition_dwell_time_us_);disp(str_msg);
str_msg=sprintf('receiver_noise_bandwidth_ %f   ', receiver_noise_bandwidth_);disp(str_msg);
str_msg=sprintf('noise_bw_scale_factor_ %f   ', noise_bw_scale_factor_);disp(str_msg);



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


triangular_covariance_matrix_bw=triangular_covariance_matrix*noise_bw_scale_factor_;


%% transfert des donnees dans une nouvelle matrices


if (oversampling==1)    
    data_before_pre_whitening= zeros(enc_Nx*2, length(meas.data),  nCoils);
else
    data_before_pre_whitening= zeros(enc_Nx, length(meas.data),  nCoils);
end

size(data_before_pre_whitening)

for p = 1:length(meas.data)
    
    data_before_pre_whitening(:,p,:)=double(meas.data{p});
    
end

%%  pre-whitening

% code issu de la fonction

in=data_before_pre_whitening;
orig_size = size(in);
nelements = prod(orig_size)/nCoils;

in = permute(reshape(in,nelements,nCoils),[2 1]);
out = reshape(permute(triangular_covariance_matrix_bw*in,[2 1]),orig_size);

end

