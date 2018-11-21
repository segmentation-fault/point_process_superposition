%   Copyright (C) 2017  Antonio Franco
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [E,EMC,EMCd] = test_superposition(genFun, CDFfun, n, k, nSamples,nTries)
% function [E,EMC,EMCd] = test_superposition(genFun, CDFfun, n, k, nSamples,nTries)
% returns the k-th moment E of the superposition of n i.i.d. point processes
% with interval distribution described by the CDF CDFfun using the Cox
% method, versus the Montecarlo superposition, done with nSamples samples
% nTries times, which returns the average k-th moment EMC along with its
% confidence interval (95%) EMCd, with samples generated with the function
% genFun. CDFfun is a function of one variable. genFun is a
% function that accepts a number of samples and generate them.

[EMC, EMCd] = mc_superimposed(genFun,n,k,nSamples,nTries);
E = num_superimposed(CDFfun,n,k);

end