% MATLAB File : difference_derivative_periodic.m
% [d] = difference_derivative_periodic(x,y,k,interval,varargin)
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Fri 12 Jun 2009 03:17:01 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Computes the k'th order finite difference derivative approximation
%   on the unstructured nodal inputs (x,y). Because periodicity is assumed here,
%   we need the input interval, a 2-vector specifying the interval of
%   approximation. The output satisfies size(d) = size(x). 
%
%   The optional input r is fed right into difference_stencil and has the same
%   meaning as in that function. 

function[d] = difference_derivative_periodic(x,y,k,interval,varargin)

global common;
prevpath = addpaths(common.bases.d1.newton.base);

xmin = interval(1); xmax = interval(2);

% Create stencil
n = length(x);
if isempty(varargin)
  [stencil,stencil_periodicity] = ...
                  difference_stencil(n,k,[],true);
else
  [stencil,stencil_periodicity] = ...
                  difference_stencil(n,k,varargin{1},true);
end

% Compute x values
XInput = x(stencil);
inds = stencil_periodicity==1;
% For indices that wrap down to 1:
XInput(inds) = xmax + (XInput(inds) - xmin);

inds = stencil_periodicity==-1;
% For indices that wrap up to n:
XInput(inds) = xmin - (xmax - XInput(inds));

% Use stencil to compute interpolants
dd = divided_difference(XInput.',y(stencil.'));

% Differentiate and evaluate the interpolants
d = newton_derivative_evaluate(XInput.',dd).';

path(prevpath);
