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

%Tests the second moment of various i.i.d. interval distributions
%superimposed n times, where n varies, cox vs Montecarlo

clc
clear
close all
reset(symengine)

tic

n = 2:10;
k = 2;
nSamples = 1e4;
nTries = 5;

q = 1;

%Gamma
shap = 5;
scal = 0.5;
D{q}.distName = ['Gamma shape = ' num2str(shap) ' scale = ' num2str(shap)];
D{q}.genFun = @(x) gamrnd(shap,scal,1,x);
D{q}.CDFfun = @(x) gamcdf(x,shap,scal);
q = q + 1;

%Exp
rate = 1;
D{q}.distName = ['Exp rate = ' num2str(rate)];
D{q}.genFun = @(x) exprnd(1/rate,1,x);
D{q}.CDFfun = @(x) expcdf(x,1/rate);
q = q + 1;

%Rayleigh
rate1 = 1;
D{q}.distName = ['Rayleigh scale = ' num2str(1/rate1)];
D{q}.genFun = @(x) raylrnd(1/rate1,1,x);
D{q}.CDFfun = @(x) raylcdf(x,1/rate1);
q = q + 1;

%Chi square
v1 = 4;
D{q}.distName = ['Chi square v = ' num2str(v1)];
D{q}.genFun = @(x) chi2rnd(v1,1,x);
D{q}.CDFfun = @(x) chi2cdf(x,v1);
q = q + 1;

%Weibull
a1 = 0.15;
b1 = 0.8;
D{q}.distName = ['Weibull a = ' num2str(a1) ' b = ' num2str(b1)];
D{q}.genFun = @(x) wblrnd(a1,b1,1,x);
D{q}.CDFfun = @(x) wblcdf(x,a1,b1);
q = q + 1;

for q = 1:numel(D)
    V = zeros(1,numel(n));
    VMC = zeros(1,numel(n));
    VMCd = zeros(1,numel(n));
    
    d0 = D{q};
    
    parfor i = 1:numel(n)
        [E,EMC,EMCd] = test_superposition( d0.genFun, d0.CDFfun, n(i), k, nSamples,nTries);
        V(i) = E;
        VMC(i) = EMC;
        VMCd(i) = EMCd;
    end
    
    if usejava('jvm') && ~feature('ShowFigureWindows')
        disp('Done');
    else
        figure;
        hold on
        plot(n,V,'x--','DisplayName','Analytical','LineWidth',2);
        errorbar(n,VMC,VMCd,'+:','DisplayName','MC','LineWidth',2);
        legend('show');
        title(D{q}.distName);
        xlabel('n');
        ylabel(['E(D^' num2str(k) ')']);
        ylim([0 inf]);
        grid on
        hold off
    end
end

elapsedTime = toc;
s = seconds(elapsedTime);
s.Format = 'hh:mm:ss.SSS';
disp('Elapsed time (hms): ');
disp(s);