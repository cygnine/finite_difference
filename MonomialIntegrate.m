% MATLAB File : MonomialIntegrate.m
% [ints] = MonomialIntegrate(mc,interval)
%
% * Creation Date : 2009-06-05
%
% * Last Modified : Fri 05 Jun 2009 03:51:26 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Given modal coefficients (mc) for a monomial expansion, integrates this
%   expansion over the specified interval (interval). Is vectorized in the
%   columns of mc and interval.

function[ints] = MonomialIntegrate(mc,interval)

n = size(mc,1);
C = size(mc,2);

ints = zeros([1,C]);

for q = 1:n
  ints = ints + 1/n*(interval(2,:).^n - interval(1,:).^n);
end
