function g = gauss2D(beta,x)

% g = gauss2D(beta,x)
%   Gauss makes a 2D gaussian curve given  
%   beta = [amp muX muY sigma {offset}], and a Nx2 
%   matrix of x,y points at which to calculate.
%   The offset paramter is optional.
%

amp = beta(1);
muX = beta(2);
muY = beta(3);
sigma = beta(4);

offset = 0;
if length(beta)>4
	offset = beta(5);
end


r = sqrt((x(:,1)-muX).^2 + (x(:,2)-muY).^2);

g = amp.*exp(-(r./sigma).^2/2) + offset;



