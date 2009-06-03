% MATLAB File : DerivativeMatrix.m
% DerivativeMatrix(x,k,{r},{periodic})
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Wed 03 Jun 2009 04:47:38 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Creates a sparse finite-difference matrix of order k on the 1D
%   mesh defined by the nodal locations x. The optional inputs r and periodic
%   are the shift and periodicity flags, respectively, and serve the same
%   purpose as in DifferenceStencil, where they are explained. 

function[mat] = DerivativeMatrix(x,k,varargin);

% Create stencil
n = length(x);
stencil = DifferenceStencil(n,k,varargin{:});

mat = spalloc(n,k+1,n*(k+1));

% This is the really slow, trivial way to do things
for q = 1:n
  y = zeros([n,1]);
  y(q) = 1;

  % Use stencil to compute interpolants
  dd = DividedDifference(x(stencil).',y(stencil.'));

  % Differentiate and evaluate the interpolants
  mat(:,q) = sparse(NewtonDiffEval(x(stencil).',dd).');
end
