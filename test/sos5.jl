using DynamicPolynomials, MomentTools

X = @polyvar x y

f = x + y + 3
G = [ y ]
H = [ x^2-1, y^2-x-2 ]

d = 2

L  = monomials(X, 0:d)

optimizer = MMT[:optimizer]

WS0, WS, P, v, M = sos_decompose(f, G, H, X, d; exact=true, round = 1 )

eq0 = f - wsos(WS0) - dot(G, wsos.(WS)) - dot(H,P)
println("reminder === ", eq0 )

Q0 = wsos(WS0)
Q1 = wsos(WS[1])

WS0, Q0
