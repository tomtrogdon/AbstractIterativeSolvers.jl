using AbstractIterativeSolvers, LinearAlgebra
using Test

@testset "AbstractIterativeSolvers.jl: GMRES" begin
    n = 20
    A = randn(n,n)/(2*sqrt(n)) + I
    b = randn(n,1)
    out = GMRES(x -> A*x, b, (x,y) -> (y'*x)[1], 1e-10, n+1, x -> x)
    sol = sum([out[2][i]*out[1][i] for i=1:length(out[2]) ])
    @test norm(A*sol-b) < 1e-10
end
