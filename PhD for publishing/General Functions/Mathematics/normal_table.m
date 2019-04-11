function [ X_t ] = normal_table( X , IC , n )

% This functions aims at creating a table of n elements from one data
% consisting in values normally distributed centered on X with a 99%
% interval of confidence equalt [X - IC ; X + IC]

%% INPUTS

% X: value to create the table from
% IC : 99% IC additive boundary aka 3*standard deviation if X is a normally distributed statistics
% n :  size of the table to create

%% OUTPUTS

% X_t : table of values created

%% Calculations

A = randn(n,1);

X_t = X + IC/3*A;

end