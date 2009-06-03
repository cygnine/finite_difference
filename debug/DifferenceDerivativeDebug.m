% MATLAB File : DifferenceDerivativeDebug.m
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Wed 03 Jun 2009 04:49:15 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Debugging file for DifferenceDerivative.

clear
cd ..

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

%%%%%%%%%%%%%% Try derivative evaluation of known function
f = @(x) sin(x);
df = @(x) cos(x);

N = 100;
xmin = 0;
xmax = 2*pi;
x = sort((xmax-xmin)*rand([N,1]) + xmin);
fx = f(x);

k = 5;

dfx = DifferenceDerivative(x,fx,k,[],false);

fprintf('Derivative error is %3.3e\n', norm(dfx-df(x)));

% While we're here, let's test making the sparse differentiation matrix, too
spmat = DerivativeMatrix(x,k,[],false);

dfx2 = spmat*fx;

fprintf('Derivative error is %3.3e\n', norm(dfx2-df(x)));

cd debug
