function [P Xp Zpred R] = fitQuadPeak(X,Zin)

%
% [P Xp Zfit R] = fitQuadPeak(X,Z);
%
%	Fits a gaussian to a funciont Z(x,y)
% Input:
%		X - Nx2 matix of locations [x y]
% 	Z - Nx1 vector of map values
% Output:
%   P - fitted parameters in log space
%   Xp - x,y locations of the gaussian peak
%   Zfit - values of the fitted gaussian to Z at locations X
%   R - correlation of the fit (sqrt of the model skill, r^2)
%



% make weights for a better peak fit
%Xpm = mean(X); %-P(4)/(P(7)*2); % first guess peaks?
%N = numel(Z);
%Xp(2) = -P(8)/(P(9)*2);
%W = ones(size(Z(:))); %try linear weights for now %sqrt(sum((X-repmat(Xp,size(X,1),1)).^-2,2)); %distance^-2 weight %
%W = Z;
%W(isinf(W)) = max(W(~isinf(W)));

% make agregate matrix
%A = repmat(W,1,9).*[ones(size(X,1),1) X(:,2) X(:,2).^2 X(:,1) X(:,1).*X(:,2) X(:,1).*X(:,2).^2 X(:,1).^2 X(:,1).^2.*X(:,2) X(:,1).^2.*X(:,2).^2];
%A = repmat(W,1,4).*[ones(size(X,1),1) X(:,1) X(:,2) (X(:,1).^2 + X(:,2).^2)];
%P = inv(A'*A)*(A'*(Z.*W));

Z = log(Zin); % gaussian fit
A = [ones(size(X,1),1) X(:,1) X(:,2) (X(:,1).^2 + X(:,2).^2)]; % no weighting
P = (A'*A)\(A'*(Z));%inv(A'*A)*(A'*(Z));
Zpred = exp(A*P);
%Zpred = A*P;

Xp(1) = -P(2)/(2*P(4));
Xp(2) = -P(3)/(2*P(4));

%% iterate to a solution try 9 steps
%Xp = Xpm;
%for j = 1:9
%	Xp(1) = -(P(4) + P(5)*Xp(2) + P(6)*Xp(2)^2)/(2*P(7) + 2*P(8)*Xp(2) + 2*P(9)*Xp(2)^2);
%	Xp(2) = -(P(2) + P(5)*Xp(1) + P(8)*Xp(1)^2)/(2*P(3) + 2*P(6)*Xp(1) + 2*P(9)*Xp(1)^2);
%end

% model correlation
mnZ = mean(Zin);
SSt = sum((Zin - mnZ).^2);
SSd = sum((Zpred - mnZ).^2);
R = sqrt(SSd/SSt);  %  correlation

