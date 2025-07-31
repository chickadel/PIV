function [outU outV cor2 R snr xc yc] = solvePIV(target,I,winN,x,y,varargin)

%
% [outU outV cor2 R snr xc yc] = solvePIV(I1,I2,winN,x,y,{method});
%
% I1 image is at t and I2 the image at t+dt
%

% test data, a simple shift down and across
% junk = rand(120); target = junk(2:101,2:101)+0.1*rand(100); I = junk(1:100,1:100); winN = 16;
% target= imrotate(junk,-2.5,'bicubic','crop'); I = imrotate(junk,2.5,'bicubic','crop');
% target = junk(2:101,2:51); I = junk(1:100,1:50); winN = 16;


% get inputs
multiFlag = 0;
if ~isempty(varargin)
	method = varargin{1};
else
	method = 3;
end
if length(varargin)>1
	offset = varargin{2};
	xi = offset{1};
	yi = offset{2};
	Uo = offset{3};
	Vo = offset{4};
	multiFlag = 1;
end

% get sizes and analysis centers
[m n] = size(I); % Image and target should be the same size!
Nt = floor([m n]/(winN/2));
%[nC mC] = meshgrid(winN:winN/2:(Nt(2)*winN/2-winN/2),winN:winN/2:(Nt(1)*winN/2-winN/2)); % 50% overlap
[nC mC] = meshgrid(winN/2:winN/2:(Nt(2)*winN/2-winN/2),winN/2:winN/2:(Nt(1)*winN/2-winN/2)); % 50% overlap
Ntotal = numel(mC);
xc = x(mC(:,1),nC(1,:));
yc = y(mC(:,1),nC(1,:));
% make non-analysis mask
maskN = interp2(x,y,target.*I,xc,yc,'*linear',nan);

% interp offsets
if multiFlag
	if all(size(Uo) < 3)
	  Uo = round(repmat(nanmedian(Uo(:)),size(xc))); % interp2 wont work with less than 2 inputs, just use the median
	  Vo = round(repmat(nanmedian(Vo(:)),size(xc)));
		%Uo = round(interp2(xi,yi,Uo,xc,yc,'*nearest',nanmedian(Uo(:))));
		%Vo = round(interp2(xi,yi,Vo,xc,yc,'*nearest',nanmedian(Vo(:))));
	else
		Uo = round(interp2(xi,yi,Uo,xc,yc,'*linear',nanmedian(Uo(:))));  %cubic
		Vo = round(interp2(xi,yi,Vo,xc,yc,'*linear',nanmedian(Vo(:))));  %cubic
	end
	%mUo = nanmedian(Uo(:));
	%mVo = nanmedian(Vo(:));
	%Uo = colfilt(Uo-mUo,[5 5],'sliding',@nanmedian)+mUo;
	%Vo = colfilt(Vo-mVo,[5 5],'sliding',@nanmedian)+mVo;
	%Uo = Uo+round(winN/8); % add in an enforced offset to get away from 0 peak
	%Vo = Vo+round(winN/8);
else
	Uo = 0*xc;
	Vo = Uo;
end

% loop to find the u and v for each center
cor2 = nan*nC;
R = cor2;
snr = cor2;
outV = cor2; % horizontal image velocity
outU = outV;
if method == 1
	opts = optimset('display','off');%optimset('display','off','tolx',10^-2,'tolfun',10^-4);%,'maxiter',1000);
else
  opts = [];
end
%smoothing taper for phaseCorr2
smoothF = ifftshift(gausswin(winN*2-1,3)*(gausswin(winN*2-1,3)')); % get rid of noise by smoothing, been using sig = 3
%smoothF = ifftshift(gausswin(winN*2-1,4)*(gausswin(winN*2-1,4)'));
smoothF([1 end],[1 end]) = 0; % set zero f to no energy

[SJ SI] = meshgrid(-(-(winN-1):(winN-1)),-(-(winN-1):(winN-1)));
itDsip = ceil(Ntotal/10);

for j = find(~isnan(maskN(:)))' %1:Ntotal %parfor
%try
	if ~mod(j,itDsip) % zero out for parfor
		fprintf('  %d of %d\r',j,Ntotal)
	end
	% set search window. make it two times as large as the target on a side, for now
	%sJ = [(nC(j)-winN+1):(nC(j)+winN)]; %[1:n]; %
	%sI = [(mC(j)-winN+1):(mC(j)+winN)]; %[1:m]; %

	% set target/template window
	tI = [mC(j)-winN/2+1:mC(j)+winN/2];
	tJ = [nC(j)-winN/2+1:nC(j)+winN/2];

	% do correlation and get rid of edge effects
	%%Ct = normxcorr2(target(tI,tJ)-mean2(target(tI,tJ)),I(sI,sJ)-mean2(I(sI,sJ)));
	%Ct = normxcorr2((target(tI,tJ)-mean2(target(tI,tJ))).*tTaper,(I(sI,sJ)-mean2(I(sI,sJ))).*sTaper);
	%Ct = Ct(winN:end-winN+1,winN:end-winN+1);
	%[Ct snr(j)] = phaseCorr2((I(tI,tJ)-mean2(I(tI,tJ)))+1,(target(tI,tJ)-mean2(target(tI,tJ)))+.1);

	if any((tI-Vo(j))<=0) | any((tJ-Uo(j))<=0) | any((tI-Vo(j))>m) | any((tJ-Uo(j))>n) % window is outside of the analysis region
	  continue
  end

	[Ct snr(j)] = phaseCorr2(I(tI,tJ),target(tI-Vo(j),tJ-Uo(j)),0,smoothF);
	%[Ct snr(j)] = phaseCorr2(I(tI+mC(j),tJ+nC(j)),target((tI+mC(j))-Vo(j),(tJ+nC(j))-Uo(j)),0,smoothF);
	Ct(winN,winN) = nan; %0;
	% find the peak location
	%[SJ SI] = meshgrid(-winN/2:winN/2,-winN/2:winN/2);
	[cor2(j) maxInd] = nanmax(Ct(:));

	switch method
	case 1
		% nonlinear gaussian fit to the correlation surface
		mask = abs(SJ-SJ(maxInd))<3 & abs(SI-SI(maxInd))<3; %get pixels around the peak
		beta0 = [cor2(j) SJ(maxInd) SI(maxInd) 3 0];
		ub = [1 beta0(2:3)+winN/4 winN/4 inf];
		lb = [0 beta0(2:3)-winN/4 0 0];
		beta = lsqcurvefit('gauss2D',beta0,[SJ(mask(:)) SI(mask(:))],Ct(mask(:)),lb,ub,opts); %[beta,rnorm,res,eflag,outp,lambda,jacb]
		%ci = diff(nlparci(beta,res,'jacobian',jacb)'); %
		outU(j) = beta(2);
		outV(j) = beta(3);
		outUc(j) = beta(4);
		outVc(j) = beta(4);
		cor2(j) = beta(1);
	case 2 % straight peak
		outU(j) = SJ(maxInd);
		outV(j) = SI(maxInd);
	case 3 % quad fit corr peak
		Ct = Ct-median(Ct(~isnan(Ct(:))));	% nanmedian
		%mask = abs(SJ-SJ(maxInd))<winN/4 & abs(SI-SI(maxInd))<winN/4 & abs(SJ)<winN/2 & abs(SI)<winN/2; %get pixels around the peak and not too far from the center
		mask = abs(SJ-SJ(maxInd))<4 & abs(SI-SI(maxInd))<4 & Ct>0; % <5 for real data % & ((SJ-SJ(maxInd)) & (SI-SI(maxInd))) %winN/4 <2
		[P Xp Cpred R(j)] = fitQuadPeak([SJ(mask(:)) SI(mask(:))],Ct(mask(:))); %fit to the peak
		%[P Xp Cpred R(j)] = fitGaussPeak([SJ(mask(:)) SI(mask(:))],Ct(mask(:))); %fit to the peak
		%outU(j) = nansum(SJ(mask(:)).*Ct(mask(:)));
		%outV(j) = nansum(SI(mask(:)).*Ct(mask(:)));
		outU(j) = Xp(1);
		outV(j) = Xp(2);
		%cor2(j) = max(Cpred(:));
	end

	% diagnostic plotting, switched off by default
	if 0
	clf
	plot(Ct,'b')
	hold on
	junk = nan*Ct;
	junk(mask) = Cpred;
	plot(junk,'r.-')
	drawnow
	end

	if 0 % comment in for 2D visualizing
	figure(2)
	clf
		subplot(121)
			imagesc([max(SJ(:)) min(SJ(:))],[max(SI(:)) min(SI(:))],Ct)
			hold on
			plot(outU(j),outV(j),'ko','markersi',14,'linew',2)
			plot(SJ(maxInd),SI(maxInd),'k+','markersi',14,'linew',1)
			axis image
	end
	if 0 % comment in for 2D visualizing
	subplot(122)
		junk = nan*Ct;
		junk(mask) = Cpred;
			imagesc([max(SJ(:)) min(SJ(:))],[max(SI(:)) min(SI(:))],junk)
			%imagesc([max(SJ(:)) min(SJ(:))],[max(SI(:)) min(SI(:))].*mask)
			hold on
			plot(outU(j),outV(j),'ko','markersi',14,'linew',2)
			plot(SJ(maxInd),SI(maxInd),'k+','markersi',14,'linew',2)
			axis image
		drawnow
	end
	if 0
	subplot(122)
			Cg = gauss2D(beta,[SJ(:) SI(:)]);
			imagesc([max(SJ(:)) min(SJ(:))],[max(SI(:)) min(SI(:))],reshape(Cg,size(Ct)).*mask)
			hold on
			plot(beta(2),beta(3),'ko','markersi',14,'linew',2)
			plot(SJ(maxInd),SI(maxInd),'k+','markersi',14,'linew',2)
			axis image
		drawnow
	end

%catch
%disp(lasterr)
%end
end

% account for multipass offset
if multiFlag
	outU = outU+Uo;
	outV = outV+Vo;
end


