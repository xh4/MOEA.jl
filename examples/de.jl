using MOEA

method = DE()

result = optimize(sphere,
                  method,
		  population = [Individual([2.0, 2.0]) for i in 1:method.populationSize],
		  options = MOEA.Options(store_trace=true))
