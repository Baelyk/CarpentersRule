function chordarcconstant = chordarc(polygon)
	n = size(polygon, 1);
	s = 1 : n - 1;
	t = 2 : n;
	d = vecnorm((polygon(s, :) - polygon(t, :))');
	G = graph(s, t, d);

	s = 1;
	a = 5;
	b = 6;

	x = polygon(a, :)';
	y = polygon(b, :)';
	z = polygon(s, :)';
	[~, lx] = shortestpath(G, s, a)
	[~, ly] = shortestpath(G, s, b)
	[~, ledge] = shortestpath(G, a, b)

	%A = norm(y - x)^2
	%B = -4 * x' * (y - x) - 2 * norm(y - x)^2 * (lx + 1)
	%C = -2 * x' * (y - x) * (lx + 1) - norm(x - s)^2

	%sols = (2 * A)^-1 * (-B + [-1, 1] * sqrt(B^2 - 4 * A * C))

	%t = (x + sols .* (y - x))'


	lambda = sym("lambda", "real")
	t = x + lambda * (y - x)
	slx = sym(lx)
	sly = sym(ly)
	sledge = sym(ledge)

	% Arc distance from z to t is min(sltx, slty), because connecting one of t -> x -> z (sltx) or t -> y -> z (slty) is shorter
	sltx = sledge * lambda + slx;
	slty = sledge * (1 - lambda) + sly;

	euclid = norm(t - z)

	f_sltx = sltx / euclid
	df_sltx = diff(f_sltx)
	f_slty = slty / euclid
	df_slty = diff(f_slty)
	lam_x = solve(df_sltx == 0, lambda)
	lam_y = solve(df_slty == 0, lambda)
	double(lam_x)
	double(lam_y)

	x
	y
	t = double(subs(t, [lam_x, lam_y])')

	clf;
	hold on;
	plot(polygon(:, 1), polygon(:, 2), "-ob");
	plot(t(:, 1), t(:, 2), "or");
end
