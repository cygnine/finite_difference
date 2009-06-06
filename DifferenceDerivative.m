% MATLAB File : DifferenceDerivative.m
% [d] = DifferenceDerivative(x,y,k,{r=0})
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Fri 05 Jun 2009 09:47:32 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Computes the k'th order finite difference derivative approximation
%   on the unstructured nodal inputs (x,y). The output satisfies size(d) =
%   size(x).
%
%   The two optional inputs r and periodic are fed right into DifferenceStencil
%   and have the same meaning as in that function. 

function[d] = DifferenceDerivative(x,y,k,varargin)

global common;
prevpath = addpaths(common.bases.d1.newton.base);

% Create stencil
n = length(x);
[stencil] = DifferenceStencil(n,k,varargin{:});

% Use stencil to compute interpolants
dd = DividedDifference(x(stencil).',y(stencil.'));

% Differentiate and evaluate the interpolants
d = NewtonDiffEval(x(stencil).',dd).';

path(prevpath);
