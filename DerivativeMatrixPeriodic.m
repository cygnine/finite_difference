% MATLAB File : DerivativeMatrixPeriodic.m
% [mat] = DerivativeMatrixPeriodic(x,k,interval,{r})
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Sat 06 Jun 2009 03:17:52 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Creates a sparse finite-difference matrix of order k on the 1D
%   mesh defined by the nodal locations x with periodic continuation over the
%   interval specified by interval. The optional input r is the shift and serves
%   the same purpose as in DifferenceStencil, where it is explained. 

function[mat] = DerivativeMatrixPeriodic(x,k,interval,varargin);

global common;
prevpath = addpaths(common.bases.d1.newton.base);

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

% Allocation
mat = spalloc(n,n,n*(k+1));

% This is the really slow, trivial way to do things
for q = 1:n
  y = zeros([n,1]);
  y(q) = 1;

  % Use stencil to compute interpolants
  dd = DividedDifference(XInput.',y(stencil.'));

  % Differentiate and evaluate the interpolants
  mat(:,q) = sparse(NewtonDiffEval(XInput.',dd).');
end

path(prevpath);
