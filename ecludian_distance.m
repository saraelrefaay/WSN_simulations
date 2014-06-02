function [ distance ] = ecludian_distance( input1X, input2X, input1Y, input2Y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    distance = sqrt((input1X - input2X)^2 +(input1Y - input2Y)^2 );
end

