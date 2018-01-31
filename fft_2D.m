function [ out ] = fft_2D( in )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

size(in)

% Reconstruct in x 
K = fftshift(fft(fftshift(in,1),[],1),1);

out=K;

% Reconstruct in y 
out = fftshift(fft(fftshift(out,2),[],2),2);



end

