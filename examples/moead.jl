using MOEA

zdt, bounds, truepf = ZDT1()

dtlz, bounds, truepf = DTLZ1()

cts = BoxConstraints(bounds[1,:], bounds[2,:])

method = NSGAII(populationSize=100)
method = NSGAII(populationSize = 100)
method = SMSEMOA(N=100)
method = MOEADAWA(N=100)

result = optimize(dtlz,
                  method,
                  constraints = cts,
		  population = [Individual([rand() for d in 1:12]) for i in 1:100],
                  options = MOEA.Options(iterations=100, store_trace=true))


method = MOEAD(N = 100)

result = optimize(zdt,
                  method,
                  constraints = cts,
		  population = [Individual([rand() for d in 1:30]) for i in 1:100],
                  options = MOEA.Options(iterations=100, store_trace=true))

