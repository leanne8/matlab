function [rts, info] = cubic(C)
tol = 10^-6;
max = 10^6;
i = 0;
if(length(C) == 0)
	info = sprintf('no polynomial is given');

elseif (length(C) == 1)
	info = sprintf('no rts');

elseif (length(C) == 2)
	a = C(1);
	b = C(2);
	firstr = -b/a;
	info = sprintf('the only rt is ', firstr);

% perform the quadratic formula
elseif (length(C) == 3)
	a = C(1);
	b = C(2);
	c = C(3);
	%firstr = (-b + sqrt(b^2 - 4*a*c))/(2*a);
    %secondr = (-b - sqrt(b^2 - 4*a*c))/(2*a);
    firstr = (2*c) / (-b - sqrt(b^2 - 4*a*c));
    secondr = ( - b - sqrt(b^2 - 4*a*c)) / (2*a);
    info = sprintf('first rt is ', firstr, ' and the second rt is ', secondr);

% multiple roots
else
	a = C(1);
	b = C(2);
	c = C(3);
	d = C(4);
	% pick a random point
	p0 = randi([-10, 10]);
	% within the 100 iterations, perform Newtons method
	while(i < 100)
		p = p0 - (polyval(C, p0)/polyval(polyder(C), p0)
		if abs(p - p0) < tol
			fprintf(p);
			break
		end
		i = i + 1
		p0 = p
	end
	% apply bisection method
	% initial guess for a and b
	a = randi([0, 100]);
	b = randi([-100, -1]);
	while i < max;
		interval = polyval(C, a) * polyval(C, b);
		avg = (a+b)/2;
		if (polyval(C, avg) == 0)
			sprintf('rt is ', avg);
			break
		end

		if interval > 0
			a = avg;
		else 
			b = avg;
		end 
	end
	




end