% MATLAB File : derivative_matrix.m
% [mat] = derivative_matrix(x,k,varargin);
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Fri 12 Jun 2009 03:12:17 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Creates a sparse finite-difference matrix of order k on the 1D
%   mesh defined by the nodal locations x. The optional inputs r and periodic
%   are the shift and periodicity flags, respectively, and serve the same
%   purpose as in difference_stencil, where they are explained. 

function[mat] = derivative_matrix(x,k,varargin);

global common;
prevpath = addpaths(common.bases.d1.newton.base);

% Create stencil
n = length(x);
stencil = difference_stencil(n,k,varargin{:});

mat = spalloc(n,k+1,n*(k+1));

% This is the really slow, trivial way to do things
for q = 1:n
  y = zeros([n,1]);
  y(q) = 1;

  % Use stencil to compute interpolants
  dd = divided_difference(x(stencil).',y(stencil.'));

  % Differentiate and evaluate the interpolants
  mat(:,q) = sparse(newton_derivative_eval(x(stencil).',dd).');
end

path(prevpath);
