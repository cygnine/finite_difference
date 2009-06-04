% MATLAB File : DifferenceStencil.m
% [stencil,{StencilPeriodicity}] = DifferenceStencil(n,k,{r=0},{periodic=false})
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Thu 04 Jun 2009 06:00:08 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Creates indices corresponding to the `default' finite-difference
%   stencil. The order k must be greater than or equal to 0 and indicates the
%   polynomial order of approximation. The output is an n x (k+1) stencil matrix
%   containing indices of an n-vector over which the difference approximation
%   is taken. Assumes that the external n-vector is ordered.
%   
%   If no additional inputs are given, this uses a central difference stencil
%   with a preference to the left in the case of odd k.
% 
%   The first optional input r (length n) indicates the index offset from the
%   central stencil. r cannot be so large that it removes the point of
%   differentiation from the stencil. An error is not raised in this case, but
%   the resulting stencil is moved so that it intersects with the point of
%   differentiation.
%  
%   The second optional input periodic is a boolean indicating whether the input is
%   periodic. If it is not, then one-sided stencils are used near the
%   boundaries.  (Default is false)
%   If the periodicity flag is set to true, then the output StencilPeriodicity
%   is a size(stencil) int8 array with values 0, \pm 1. +1 indicates that the
%   nodal index had a value greater than n and was wrapped down, and -1
%   indicates that the nodal index had a value less than 1 and was wrapped up. 

function[stencil,varargout] = DifferenceStencil(n,k,varargin)

% Input data parsing
if length(varargin)==2
  periodic = varargin{2};
  r = varargin{1};
  if isempty(r)
    r = zeros([n,1]);
  end
elseif length(varargin)==1
  r = varargin{1};
  periodic = false;
else
  r = zeros([n,1],'int8');
  periodic = false;
end
if length(r)==1
  r = r*ones([n,1],'int8');
end
r = int8(r);

% First let's just create the linear offsets, ignoring boundaries
stencil = zeros([n,k+1],'int8');

% Lots of initial data is needed
stencil(:,1) = 1:n;
NegativeK = true([n,1]);
NegativeCount = zeros([n,1],'int8');
PositiveCount = zeros([n,1],'int8');
Rsign = sign(r);
Rmag = abs(r);

% Downgrade r to its maximum values
if mod(k,2)==0
  Rmag = min(Rmag,k/2);
else
  Rmag = min(Rmag,(k+1)/2);
  inds = and(Rmag>(k-1)/2, Rsign<0);
  Rmag(inds) = (k-1)/2;
end
r = Rmag.*Rsign;

% There's probably a better way to do this
for q = 1:k
  % Do the default left-stencil creation:
  inds = NegativeK;
  stencil(inds,q+1) = stencil(inds,1) -NegativeCount(inds) -1;
  NegativeCount(inds) = NegativeCount(inds) + 1;

  % Do the default right-stencil creation:
  inds = not(NegativeK);
  stencil(inds,q+1) = stencil(inds,1) + PositiveCount(inds)+1;
  PositiveCount(inds) = PositiveCount(inds) + 1;

  % Update NegativeK 
  NegativeK = not(NegativeK);
end
% Shift by r
stencil(:,2:(k+1)) = stencil(:,2:(k+1)) + repmat(r,[1,k]);
% Find which rows have repeated indices
inds = stencil(:,2:(k+1))==repmat(stencil(:,1),[1,k]);

stencil(:,2:(k+1)) = stencil(:,2:(k+1)) + ...
         int8(inds).*repmat(r,[1,k]);

% Deal with boundaries: again, very stupid, but it works
if not(periodic)
  % Right boundary:
  inds = stencil>n;
  rows = find(any(inds,2));
  for q = rows'
    RowInds = inds(q,:);
    cols = sort(find(RowInds));
    MinN = min(stencil(q,:));
    ColCount = 1;
    for k = cols
      stencil(q,k) = MinN - ColCount;
      ColCount = ColCount + 1;
    end
  end

  % Left Boundary
  inds = stencil<1;
  rows = find(any(inds,2));
  for q = rows'
    RowInds = inds(q,:);
    cols = sort(find(RowInds));
    MaxN = max(stencil(q,:));
    ColCount = 1;
    for k = cols
      stencil(q,k) = MaxN + ColCount;
      ColCount = ColCount + 1;
    end
  end

else  % Periodic case is *much* easier
  StencilPeriodicity = zeros([n,k+1],'int8');
  StencilPeriodicity(stencil>n) = +1;
  StencilPeriodicity(stencil<1) = -1;
  varargout{1} = StencilPeriodicity;

  stencil = mod(stencil-1,n)+1;
end
