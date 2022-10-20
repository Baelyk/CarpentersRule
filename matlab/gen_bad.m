exitflag = 0
while exitflag ~= -2
	P = polygon_fxn(4)
	[v exitflag] = find_velocity(P(1 : end - 1, :), false, 10^-3)
end
