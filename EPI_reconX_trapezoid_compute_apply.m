function [ Mpos_ ,   Mneg_] = EPI_reconX_trapezoid_compute_apply( reconx , trajectoryPos_ , trajectoryNeg_ , hdr_in)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Compute the reconstruction operator
Km = floor(reconx.encodeNx_ / 2.0);
Ne = 2*Km + 1;

% resize the reconstruction operator
Mpos_=zeros(reconx.reconNx_,reconx.numSamples_);
Mneg_=zeros(reconx.reconNx_,reconx.numSamples_);

% evenly spaced k-space locations
keven = linspace(-Km, Km, Ne);
%keven.print("keven =");

% image domain locations [-0.5,...,0.5)
x = linspace(-0.5,(reconx.reconNx_-1.)/(2.*reconx.reconNx_),reconx.reconNx_);
%x.print("x =");

% DFT operator
% Going from k space to image space, we use the IFFT sign convention
F=zeros(reconx.reconNx_, Ne);
fftscale = 1.0 / sqrt(Ne);

for p=1:1:reconx.reconNx_
    for q=1:1:Ne
        F(p,q) = fftscale * exp(complex(0.0,1.0*2*pi*keven(q)*x(p)));
    end
end
%F.print("F =");

% forward operators
Qp=zeros(reconx.numSamples_, Ne);
Qn=zeros(reconx.numSamples_, Ne);

for p=1:1:reconx.numSamples_
    %GDEBUG_STREAM(trajectoryPos_(p) << "    " << trajectoryNeg_(p) << std::endl);
    for q=1:1:Ne
        Qp(p,q) = sinc(trajectoryPos_(p)-keven(q));
        Qn(p,q) = sinc(trajectoryNeg_(p)-keven(q));
    end
end

%Qp.print("Qp =");
%Qn.print("Qn =");

% recon operators
Mp=zeros(reconx.reconNx_,reconx.numSamples_);
Mn=zeros(reconx.reconNx_,reconx.numSamples_);
Mp = F * pinv(Qp);
Mn = F * pinv(Qn);

% Compute the off-center correction:     /////

% Compute the off-center distance in the RO direction:
roOffCenterDistance = calcOffCenterDistance( hdr_in );
%GDEBUG_STREAM("roOffCenterDistance_: " << roOffCenterDistance << ";       encodeFOV_: " << encodeFOV_);

my_keven = linspace(0, reconx.numSamples_ -1, reconx.numSamples_);
%         % find the offset:
%         % PV: maybe find not just exactly 0, but a very small number?
trajectoryPosArma = trajectoryPos_;

n = find( trajectoryPosArma==0, 1, 'first');

my_keven = my_keven - (n-1);  %%attention 128

% my_keven
%         % Scale it:
%         % We have to find the maximum k-trajectory (absolute) increment:
Delta_k = abs( trajectoryPosArma(2:reconx.numSamples_,1) - trajectoryPosArma(1:reconx.numSamples_-1,1));
% max(Delta_k(:))

my_keven = my_keven *max(Delta_k(:));
%
%         % off-center corrections:


clear offCenterCorrN offCenterCorrP

myExponent = zeros(reconx.numSamples_,1);


for l=1:size(trajectoryPosArma,1)
    myExponent(l,1)=imag( 2*pi*roOffCenterDistance/reconx.encodeFOV_*(trajectoryPosArma(l,1)-my_keven(1,l)) );
end

offCenterCorrN = exp( myExponent );

for l=1:size(trajectoryPosArma,1)
    myExponent(l,1)=imag( 2*pi*roOffCenterDistance/reconx.encodeFOV_*(trajectoryNeg_(l,1)+my_keven(1,l)) );
end

offCenterCorrP = exp( myExponent );


%         %    for (q=0; q<numSamples_; q++) {
%         %      GDEBUG_STREAM("keven(" << q << "): " << my_keven(q) << ";       trajectoryPosArma(" << q << "): " << trajectoryPosArma(q) );
%         %      GDEBUG_STREAM("offCenterCorrP(" << q << "):" << offCenterCorrP(q) );
%         %    }


%         Finally, combine the off-center correction with the recon operator:
Mp = Mp * diag(offCenterCorrP);
Mn = Mn * diag(offCenterCorrN);

%          and save it into the NDArray members:
for p=1:1:reconx.reconNx_
    for q=1:1:reconx.numSamples_
        Mpos_(p,q) = Mp(p,q);
        Mneg_(p,q) = Mn(p,q);
    end
end

% size(Mpos_)
% size(Mneg_)
%Mp.print("Mp =");
%Mn.print("Mn =");


return


function roOffCenterDistance=calcOffCenterDistance( hdr_in)

% armadillo vectors with the position and readout direction:
pos=zeros(3,1);
RO_dir=zeros(3,1);

for i=1:1:3
    pos(i,1)    = hdr_in.position(i);
    RO_dir(i,1) = hdr_in.read_dir(i);
end

roOffCenterDistance = dot(pos, RO_dir);
%GDEBUG_STREAM("roOffCenterDistance: " << roOffCenterDistance );

roOffCenterDistance


return
