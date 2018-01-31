function [ trajectoryPos_ , trajectoryNeg_ , reconx , MB_factor, Blipped_CAIPI, MB_Slice_Inc] = EPI_reconX_trapezoid_compute_trajectory( reconx , temp_struct_long, temp_struct_double )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for ll=1:size(temp_struct_long,2)
    
    if (strcmp(temp_struct_long(ll).name,'numberOfNavigators'))
        numberOfNavigators_=temp_struct_long(ll).value;
        str_msg=sprintf('%s %f', temp_struct_long(ll).name , temp_struct_long(ll).value); disp(str_msg);
    end
    if (strcmp(temp_struct_long(ll).name,'rampUpTime'))
        rampUpTime_=temp_struct_long(ll).value;
        str_msg=sprintf('%s %f', temp_struct_long(ll).name , temp_struct_long(ll).value); disp(str_msg);
    end
    
    if (strcmp(temp_struct_long(ll).name,'rampDownTime'))
        rampDownTime_=temp_struct_long(ll).value;
        str_msg=sprintf('%s %f', temp_struct_long(ll).name , temp_struct_long(ll).value); disp(str_msg);
    end
    
    if (strcmp(temp_struct_long(ll).name,'acqDelayTime'))
        acqDelayTime_=temp_struct_long(ll).value;
        str_msg=sprintf('%s %f', temp_struct_long(ll).name , temp_struct_long(ll).value); disp(str_msg);
    end
    
    if (strcmp(temp_struct_long(ll).name,'numSamples'))
        numSamples_=temp_struct_long(ll).value;
        str_msg=sprintf('%s %f', temp_struct_long(ll).name , temp_struct_long(ll).value); disp(str_msg);
    end
    
    if (strcmp(temp_struct_long(ll).name,'flatTopTime'))
        flatTopTime_=temp_struct_long(ll).value;
        str_msg=sprintf('%s %f', temp_struct_long(ll).name , temp_struct_long(ll).value); disp(str_msg);
    end
    
    if (strcmp(temp_struct_long(ll).name,'flatTopTime'))
        flatTopTime_=temp_struct_long(ll).value;
        str_msg=sprintf('%s %f', temp_struct_long(ll).name , temp_struct_long(ll).value); disp(str_msg);
    end    
    
end



for ll=1:size(temp_struct_double,2)
    
     if (strcmp(temp_struct_double(ll).name,'dwellTime'))
     dwellTime_=temp_struct_double(ll).value;
     str_msg=sprintf('%s %f', temp_struct_double(ll).name , temp_struct_double(ll).value); disp(str_msg);
     end
    
       if (strcmp(temp_struct_double(ll).name,'NMBSliceBands'))
     MB_factor=temp_struct_double(ll).value;
     str_msg=sprintf('%s %f', temp_struct_double(ll).name , temp_struct_double(ll).value); disp(str_msg);
       
      end     
     
      if (strcmp(temp_struct_double(ll).name,'BlipFactorSL'))
     Blipped_CAIPI=temp_struct_double(ll).value;
     str_msg=sprintf('%s %f', temp_struct_double(ll).name , temp_struct_double(ll).value); disp(str_msg);
     
      end     
      
      if (strcmp(temp_struct_double(ll).name,'MultiBandSliceInc'))
     MB_Slice_Inc=temp_struct_double(ll).value;
     str_msg=sprintf('%s %f', temp_struct_double(ll).name , temp_struct_double(ll).value); disp(str_msg);
     
      end     
     
     
end

reconx.numSamples_=numSamples_;
reconx.rampUpTime_=rampUpTime_;
reconx.rampDownTime_=rampDownTime_;
reconx.acqDelayTime_=acqDelayTime_;
reconx.flatTopTime_=flatTopTime_;
reconx.dwellTime_=dwellTime_;

% Initialize the k-space trajectory arrays
trajectoryPos=zeros(numSamples_,1);
trajectoryNeg=zeros(numSamples_,1);

%  Temporary trajectory for a symmetric readout
%  first calculate the integral with G = 1;
nK = numSamples_;
k=zeros(nK,1);
%   float t;

%    Some timings
totTime = rampUpTime_ + flatTopTime_ + rampDownTime_;
readTime = dwellTime_ * numSamples_;

balanced_=1;
% Fix the acqDelayTime for balanced acquisitions
if (balanced_==1)
    acqDelayTime_ = 0.5 * (totTime - readTime);
    str_msg=sprintf('acqDelayTime_ %f \n', acqDelayTime_); disp(str_msg);
end

% Some Areas
totArea = 0.5*rampUpTime_ + flatTopTime_ + 0.5*rampDownTime_;
readArea =  0.5*rampUpTime_ + flatTopTime_ + 0.5*rampDownTime_;

if (rampUpTime_ > 0.0)
    readArea =  readArea - 0.5*(acqDelayTime_)*(acqDelayTime_)/rampUpTime_;
end
if (rampDownTime_ > 0.0)
    readArea = readArea - 0.5*(totTime - (acqDelayTime_+readTime))*(totTime - (acqDelayTime_+readTime))/rampDownTime_;
end

% Prephase is set so that k=0 is halfway through the readout time
prePhaseArea = 0.5 * totArea;

% The scale is set so that the read out area corresponds to the number of encoded points
scale = reconx.encodeNx_ /readArea;

for n=1:1:nK
    
    t = ((n-1)+1.0)*dwellTime_ + acqDelayTime_;  % end of the dwell time
    if (t <= rampUpTime_)
        % on the ramp up
        k(n,1) = 0.5 / rampUpTime_ * t*t;
        
    elseif ((t > rampUpTime_) && (t <= (rampUpTime_+flatTopTime_)))
        % on the flat top
        k(n,1) = 0.5*rampUpTime_ + (t - rampUpTime_);
    else
        % on the ramp down
        v = (rampUpTime_+flatTopTime_+rampDownTime_-t);
        k(n,1) = 0.5*rampUpTime_ + flatTopTime_ + 0.5*rampDownTime_ - 0.5/rampDownTime_*v*v;
    end
%     str_msg=sprintf('%d %f %f \n', n, t , k(n,1)); disp(str_msg);
end


%   // Fill the positive and negative trajectories
for n=1:1:numSamples_
    
    trajectoryPos_(n,1) = scale * (k(n,1) - prePhaseArea);
    trajectoryNeg_(n,1) = scale * (-1.0*k(n,1) + totArea - prePhaseArea);
%     str_msg=sprintf('%d %f %f \n', n,trajectoryPos_(n,1) , trajectoryNeg_(n,1)); disp(str_msg);
end


end

