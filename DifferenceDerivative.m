% MATLAB File : DifferenceDerivative.m
% DifferenceDerivative(x,y,k,{r=0},{periodic=false})
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Wed 03 Jun 2009 04:34:18 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Computes the k'th order finite difference derivative approximation
%   on the unstructured nodal inputs (x,y). 
%
%   The two optional inputs r and periodic are fed right into DifferenceStencil
%   and have the same meaning as in that function. 

function[d] = DifferenceDerivative(x,y,k,varargin)

% Create stencil
n = length(x);
stencil = DifferenceStencil(n,k,varargin{:});

% Use stencil to compute interpolants
dd = DividedDifference(x(stencil).',y(stencil.'));

% Differentiate and evaluate the interpolants
d = NewtonDiffEval(x(stencil).',dd).';
