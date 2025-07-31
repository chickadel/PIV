function [outU outV cor2 R snr xc yc] = solveMultiPassPIV(target,I,winN,x,y,varargin)

%
% [U V cor fitR SNR xc yc] = solveMultiPassPIV(I1,I2,winN,x,y,{method})
%
% Input:
%   I1,I2 - these are identically sized, equally spaced (dx and dy are equal) gidded images separated
%           by some time, dt (I1 occurs before I2)
%   winN - a vector of analysis window sizes (windows are square) for multipass PIV
%          (e.g winN = [128 64 32], would process the PIVin steps starting with a 128x128 window,
%          decreasing to 64x64, and then 32x32 points, resulting in finer resolution with each step)
%   x,y - grid of x and y locations of the data in I1 and I2, sized the same as I1 and I2
%   method - {optional} method for which peak finding algorithm is used.  Use: 1 = nonlinear search (slow),
%            2 = max phase correlatoion (fast and noisy), 3 = linear guassian fit (default).
%
% Output:
%   U, V - matricies of image horizontal (U) and vertical (V) shifts in pixerl units
%          (transform to velocities with dt and dx)
%   cor - peak phase correlation at each analysis location
%   fitR - skill of the peak estimating model fit
%   SNR - Signal to noise ratio of the analysis regions, indiciative of the available image contrast
%   xc, yc - locations of the PIV velocity grid in the coordinates system of the image (i.e. x and y)
%

try
pctRunOnAll warning off
catch
warning off
end
% get inputs
method = 3;
if ~isempty(varargin)
	method = varargin{1};
end

% run first window
fprintf(1,'                  -   pass 1 of %d, win length = %d    \r',length(winN),winN(1))
[outU outV cor2 R snr xc yc] = solvePIV(target,I,winN(1),x,y,method);


% loop through winN
for j = 2:length(winN)
	fprintf(1,'                  -   pass %d of %d, win length = %d    \r',j,length(winN),winN(j))
	%magU = abs(outU +outV*i);
	%Uthresh = nanmedian(magU(:))+ 2*nanstd(magU(:));
	%maskU = magU >= Uthresh;
	outU(isnan(outU)) = nanmean(outU(:)); % | maskU) = nanmedian(outU(:));
	outV(isnan(outV)) = nanmean(outV(:)); % | maskU) = nanmedian(outV(:));
	[outU outV cor2 R snr xc yc] = solvePIV(target,I,winN(j),x,y,method,{xc yc outU outV});
end

