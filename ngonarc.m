function points = ngonarc(n)
	% The nth roots of harmony or whatever: exp(2 * pi * i * k / n) gives kth
	points = exp(2 * pi * 1i * (0 : n - 1) / n);
	% Convert to 2D vectors
	points = [real(points)', imag(points)'];
	% Translate so all points are in the first quadran
	points = points - min(points);
end
