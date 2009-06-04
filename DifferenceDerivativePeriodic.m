% MATLAB File : DifferenceDerivativePeriodic.m
% [d] = DifferenceDerivativePeriodic(x,y,k,interval,{r=0})
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Wed 03 Jun 2009 08:03:25 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Computes the k'th order finite difference derivative approximation
%   on the unstructured nodal inputs (x,y). Because periodicity is assumed here,
%   we need the input interval, a 2-vector specifying the interval of
%   approximation. The output satisfies size(d) = size(x). 
%
%   The optional input r is fed right into DifferenceStencil
%   and has the same meaning as in that function. 

function[d] = DifferenceDerivativePeriodic(x,y,k,interval,varargin)

xmin = interval(1); xmax = interval(2);

% Create stencil
n = length(x);
if isempty(varargin)
  [stencil,StencilPeriodicity] = ...
                  DifferenceStencil(n,k,[],true);
else
  [stencil,StencilPeriodicity] = ...
                  DifferenceStencil(n,k,varargin{1},true);
end

% Compute x values
XInput = x(stencil);
inds = StencilPeriodicity==1;
% For indices that wrap down to 1:
XInput(inds) = xmax + (XInput(inds) - xmin);

inds = StencilPeriodicity==-1;
% For indices that wrap up to n:
XInput(inds) = xmin - (xmax - XInput(inds));

% Use stencil to compute interpolants
dd = DividedDifference(XInput.',y(stencil.'));

% Differentiate and evaluate the interpolants
d = NewtonDiffEval(XInput.',dd).';
