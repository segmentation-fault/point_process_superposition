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

function E = num_superimposed(fun,n,k)
%function E = num_superimposed(fun,n,k)
% returns the k-th moment E of the superposition of n i.i.d. point processes
% with interval distribution described by the CDF fun using the Cox
% method:
% On the Superposition of Renewal Processes
% Author(s): D. R. Cox and  Walter L. Smith
% Source: 
% Biometrika,
%  Vol. 41, No. 1/2 (Jun., 1954), pp. 91-99
% Published by: Oxford University Press on behalf of Biometrika Trust
% Stable URL: http://www.jstor.org/stable/2333008

syms x t positive
mu = integral(@(t) 1-fun(t),0,inf);

Fc = @(t) 1 - fun(t);

Fp = @(t) Fc(t) .* integral(@(x) Fc(x)./mu,t,inf).^(n-1);

E = double(k*integral(@(t) t.^(k-1) .* (Fp(t)),0,inf,'ArrayValued',true));

end