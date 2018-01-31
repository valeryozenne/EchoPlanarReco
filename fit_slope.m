function [ slope, intercept, x ] = fit_slope( navdata_ )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%navdata est de dimension

readout=size(navdata_,1);
% channels=size(navdata_,4);

ctemp =  zeros(readout,1);   % temp column complex
% tvec = zeros(readout,1);          % temp column real
x(:,1) = linspace(-0.5, 0.5,readout); % Evenly spaced x-space locations

temp0(:,:)=navdata_(:,1,:)+navdata_(:,3,:);
temp1(:,:)=conj(navdata_(:,1,:)+navdata_(:,3,:));
temp2(:,:)=conj(navdata_(:,1,:)+navdata_(:,3,:)).*navdata_(:,2,:) ;

ctemp=sum(temp2,2);

clear vec1 vec2

%% calculons la pente comme dans le gadgetron
vec1=ctemp(1:readout-2,1);
vec2=ctemp(2:readout-1,1);

% dot(conj(vec1),vec2)

%% difference dans le calcul
slope = (readout-1) * angle(dot(vec1,vec2));

ramp=complex(zeros( readout, 1), -slope*x( :, 1) );

ctemp_after = ctemp .* exp(ramp);

intercept = angle(sum(ctemp_after));


% close all
% 
% figure(1)
% plot(angle(ctemp(:,1)), 'r');
% hold on;
% plot(angle(ctemp_after(:,1)),'g');
% legend('avant','apres','Location','northwest');
end

