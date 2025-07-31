function [C effDOF] = phaseCorr2(A,B,varargin)

%
% [C {edof}] = phaseCorr2(A,B)
%	Make sure A and B are demeaned!
% A and B can be different sizes, but they may *not*
% have odd length dimensions (eg. *not* 3x5).

% compute SNR
mnAB = [mean(A(:)) mean(B(:))];
A = A - mnAB(1);
B = B - mnAB(2);
if nargout == 2
	%effDOF = 2*sum(abs(Fa(:))).^2./sum(abs(Fa(:)).^2);
	effDOF = sum([A(:) B(:)].^2);
	effDOF = prod(effDOF)^.25; % image std
end

% condition inputs
if numel(A)<numel(B)
	A = padarray(A,[size(B)-size(A)]/2);
end
szA = size(A);
szB = size(B);
outsize = szA + szB - 1;
smoothF = 1; % default to smoothing
if ~isempty(varargin)
	if length(varargin) == 1
		singleVal = varargin{1};
	elseif length(varargin) > 1
		singleVal = varargin{1};
		smoothF = varargin{2};
	end
else
	singleVal = 0;
end
if (numel(smoothF) == 1) && smoothF % make a smoothing window if you want it
	smoothF = ifftshift(myGaussWin(outsize(1),3)*(myGaussWin(outsize(2),3)')); % get rid of noise by smoothing
	%smoothF = ifftshift(gausswin(outsize(1),1)*(gausswin(outsize(2),1)'));
	smoothF([1:2 end-1:end],[1:2 end-1:end]) = 0; % set zero f to no energy
elseif (numel(smoothF) == 1) && ~smoothF
	smoothF = 1;  % no smoothing, but have to set gain to 1!
end

% compute fourier trans.
Fa = fftn(A,[outsize(1) outsize(2)]);
Fb = fftn(B,[outsize(1) outsize(2)]);
Fbac = Fb.*conj(Fa); % complex conjugate

% compute phase correlation;
%C = Fbac.*smoothF.*abs(Fbac)./((abs(Fa).^2).*(abs(Fb).^2)); % try coherence filter weighting, not quite right though, need to smooth!
C = Fbac.*smoothF./abs(Fbac);  % smooth at the same time!
C(isnan(C)) = 0;
C = ifftn(C,[outsize(1) outsize(2)],'symmetric');
C = fftshift(C);
C = real(C);
%C = real(fftshift(ifft2(Fb.*conj(Fa).*(smoothF)./abs(Fb.*conj(Fa)))));

% return a single corr value at 0 phase offset if called
if singleVal
	C = C(round(outsize(1)/2),round(outsize(2)/2));
end


%%%%%%%%%%%%%%%%%%%%% myGaussWin subfunction
function w = myGaussWin(N,a)
N = N-1;
n = (0:N)'-N/2;
w = exp(-(1/2)*(a*n/(N/2)).^2);


