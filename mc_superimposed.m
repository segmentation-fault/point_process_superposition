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

function [E, Ed] = mc_superimposed(fun,n,k,nSamples,nTries)
%function [E, Ed] = mc_superimposed(fun,n,k,nSamples,nTries) 
% Montecarlo superposition of n i.i.d. point processes, done with nSamples samples
% nTries times, which returns the average k-th moment EMC along with its
% confidence interval (95%) EMCd, with samples generated with the function
% fun. fun is a function that accepts a number of samples
% and generate them.

Ev = zeros(1,nTries);

for j=1:nTries
    A = zeros(n,nSamples);

    for i = 1:n
        A(i,:) = cumsum(fun(nSamples));
    end

    S = sort(A(:));

    D = diff(S);

    Ev(j) = mean(D.^k);
end

E = mean(Ev);
Ed = 1.96*std(Ev)/nTries;

end