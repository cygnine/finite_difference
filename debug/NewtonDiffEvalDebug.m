% MATLAB File : NewtonDiffEvalDebug.m
%
% * Creation Date : 2009-06-03
%
% * Last Modified : Wed 03 Jun 2009 02:28:30 PM EDT
%
% * Created By : Akil Narayan
%
% * Purpose : Debugging file for NewtonDiffEval.

clear
cd ..

RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

%%%%%%%%%%%%%% Try derivative evaluation of known polynomial
f = @(x) x.^5 + 4*x.^3;
df = @(x) 5*x.^4 + 12*x.^2;

NTrials = 8;
x = rand([6 NTrials]);
fx = f(x);

dd = DividedDifference(x,fx);
fPrime = NewtonDiffEval(x,dd);

fprintf('Derivative error is %3.3e\n', norm(fPrime-df(x(1,:))));

cd debug
