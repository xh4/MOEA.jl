@testset "DTLZ1" begin
    fn, bounds, truepf = DTLZ1()
    @test fn(fill(0.5, 7)) == [0.125, 0.125, 0.25]                              
end

@testset "DTLZ2" begin
    fn, bounds, truepf = DTLZ2()
    @test fn(fill(0.5, 12)) == [0.5, 0.5, 0.7071]                              
end

@testset "DTLZ3" begin
    fn, bounds, truepf = DTLZ3()
    @test fn(fill(0.5, 12)) == [0.5, 0.5, 0.7071]                              
end

@testset "DTLZ4" begin
    fn, bounds, truepf = DTLZ4()
                           
end

@testset "DTLZ5" begin
    fn, bounds, truepf = DTLZ5()
    @test fn(fill(0.5, 12)) == [0.5, 0.5, 0.7071]                            
end

@testset "DTLZ6" begin
    fn, bounds, truepf = DTLZ6()
    @test fn(fill(0.5, 12)) == [0.9665, 0.9665, 1.3669]                              
end