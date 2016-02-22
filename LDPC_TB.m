clear all; close all; clc;
format compact;
input = logical(randi([0 1], 324, 1));
output = LDPC(input);
H = dvbs2ldpc(1/2);
K = H(1:200,1:200);
biterr(input, output)