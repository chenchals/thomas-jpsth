function [ angle ] = convert_tgt_octant_to_angle( octant )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

NUM_TRIALS = length(octant);

angle = NaN(1,NUM_TRIALS);

angle(octant == 1) = 0;
angle(octant == 2) = pi/4;
angle(octant == 3) = pi/2;
angle(octant == 4) = 3*pi/4;
angle(octant == 5) = pi;
angle(octant == 6) = 5*pi/4;
angle(octant == 7) = 3*pi/2;
angle(octant == 8) = 7*pi/4;

end%function:convert_tgt_octant_to_angle()

