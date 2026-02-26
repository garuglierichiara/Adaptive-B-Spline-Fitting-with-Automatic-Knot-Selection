function i = findspan_periodic(p, u, n, U, a, b)

u = a + mod(u - a, b - a);   % Map u in [a,b)
% Mod ritorna il resto della divisione intera => mod(x,L) appartiene a [0,L)

% Call findspan
i = findspan(p, u, n, U);
end
