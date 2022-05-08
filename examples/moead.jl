using MOEA

zdt, bounds, pfront = ZDT1()

constraints = BoxConstraints(bounds[1,:], bounds[2,:])

method = MOEAD(N=100)

result = optimize(zdt,
                  method,
                  constraints = constraints,
		  population = [Individual([rand() for d in 1:30]) for i in 1:100],
                  options = MOEA.Options(iterations=200, store_trace=true))
