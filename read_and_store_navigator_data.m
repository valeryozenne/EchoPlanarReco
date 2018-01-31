function [ data_navigator ] = read_and_store_navigator_data( meas, acq_phasecorr , data_in)
%
% [ data_navigator ] = read_and_store_navigator_data( meas, acq_phasecorr , data_in)
%   
%   Extract navigator data and store it.
%
%   INPUT:
%       meas [kx,ky,coil]          : Source data for grappa kernel estimation (k-space)
%       acq_phasecorr [1, n]       : 
%       data_in  [kx,ky,coil]      : Source data
%
%   OUTPUT:
%       data_navigator [kx,navigatorNumber,coil]    : 


data_navigator = zeros(size(data_in,1), length(acq_phasecorr),  size(data_in,3));

for p = 1:length(acq_phasecorr)
    
    if (p==1) 
        
       mise_a_zero=acq_phasecorr(p)-1;
        
    end
    
    ky = meas.head.idx.kspace_encode_step_1(acq_phasecorr(p)) + 1;
    
    index_navigateur=acq_phasecorr(p)-mise_a_zero;
    
    str_msg=sprintf('p %d  ky %d  index navigateur %d \n',p , acq_phasecorr(p), index_navigateur); disp(str_msg);
    
    data_navigator(:,index_navigateur,:) = data_in(:,acq_phasecorr(p),:);
    
end


return

