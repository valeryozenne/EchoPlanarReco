function [meas , meas_noise ] = read_dset( dset_noise  , dset )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% TODO add the other possibilites

%% Read all the data
% Reading can be done one acquisition (or chunk) at a time,
% but this is much faster for data sets that fit into RAM.
D_noise = dset_noise.readAcquisition();
D = dset.readAcquisition();

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

nbLines_noise=size(meas_noise.data,2);

%nb de lignes
str_msg=sprintf('le nombre TOTAL de lignes dans le fichier noise est %d \n', nbLines_noise); disp(str_msg);

%% quelques vérifications



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

% pre-whitening (NoiseAdjustGadget.cpp NoiseAdjustGadget.h)

nbLines_totale=size(meas.data,2);

%nb de lignes
str_msg=sprintf('le nombre TOTAL de lignes dans le fichier data est %d \n', nbLines_totale); disp(str_msg);



end

