% MATLAB File : difference_derivative.m
% [d] = difference_derivative(x,y,k,varargin)
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Fri 12 Jun 2009 03:15:24 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Computes the k'th order finite difference derivative approximation
%   on the unstructured nodal inputs (x,y). The output satisfies size(d) =
%   size(x).
%
%   The two optional inputs r and periodic are fed right into difference_stencil
%   and have the same meaning as in that function. 

function[d] = difference_derivative(x,y,k,varargin)

global common;
prevpath = addpaths(common.bases.d1.newton.base);

% Create stencil
n = length(x);
[stencil] = difference_stencil(n,k,varargin{:});

% Use stencil to compute interpolants
dd = divided_difference(x(stencil).',y(stencil.'));

% Differentiate and evaluate the interpolants
d = newton_derivative_evaluate(x(stencil).',dd).';

path(prevpath);
