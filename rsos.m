function [ output] = rsos( input )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


 output=sqrt(sum(abs(input).^2,3));

end

