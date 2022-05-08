# MOEA

Multi-Objective Evolutionary Algorithms in Julia

![GUI](https://user-images.githubusercontent.com/2956767/167279831-2943130e-e8e3-4bba-8c41-a31089c30ea9.png)

## Algorithms

- NSGA-II
- MOEA/D
- MOEA/D-DE

## Problems

- ZDT
- DTLZ
- WFG

## Indicators

- IGD
- Hypervolume

## GA

```julia
ga = GA(populationSize=100,
		selection=susinv,
		mutation=BGA(ones(2)),
		crossover=DC)

result = optimize(sphere,
                  ga,
				  population = [Individual([2.0, 2.0]) for i in 1:ga.populationSize],
				  options = MOEA.Options(store_trace=true))
```

## MOEA/D

```julia
zdt, bounds, pareto_solutions = ZDT1()

constraints = BoxConstraints(bounds[1,:], bounds[2,:])

method = MOEAD(N=100)

result = optimize(zdt,
                  method,
                  constraints = constraints,
			      population = [Individual([rand() for d in 1:30]) for i in 1:100],
                  options = MOEA.Options(iterations=200, store_trace=true))
```


## Test

```sh
julia -e 'import Pkg; Pkg.activate(pwd()); Pkg.test()'
```

