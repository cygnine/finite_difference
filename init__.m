function[finite_difference] = init__()
% init__ -- Initialization file for finite_difference package
%
% [nodes] = init__()

finite_difference = recurse_files(pwd);
finite_difference.debug = matlab_import('debug');
