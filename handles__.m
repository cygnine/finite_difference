function[hs,pathadditions,name] = handles__()
% handles__ -- handle structure for finite_difference module
%
% [hs,pathadditions,name] = handles__()
%
%     Returns directory pointers for common module in HS. PATHADDITIONS is a
%     cell array with a string in each element indicated paths to add to the
%     global path structure. 

name = 'finite_difference';

% This is by default
hs.base = fileparts(mfilename('fullpath'));

pathadditions = cell(0);
