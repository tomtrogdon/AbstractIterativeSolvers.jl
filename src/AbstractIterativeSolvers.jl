module AbstractIterativeSolvers

using LinearAlgebra

export GMRES, GMRES_verbose

function GMRES(A,b,inner,tol,n,cond)
    nom = a -> sqrt(abs(inner(a,a)))
    H = zeros(Complex{Float64},n+1,n)
    bnorm = nom(b)
    x = 0.
    Q = [(1.0/bnorm)*b]
    for i = 1:n
       #tic()
       #println("Operator application: ")
       v = A(Q[i])
       #toc()
       #tic()
       #println("Inner products: ")
       for j = 1:i
           #tic()
           H[j,i] = inner(Q[j],v)
           #toc()
           v = cond(v - H[j,i]*Q[j])
       end
       v = cond(v)
       #println("Assembling Q:")
       H[i+1,i] = nom(v)
       Q = vcat(Q,[copy((1.0/H[i+1,i])*v)])
       #print("Arnoldi: ")
       #toc()
       #return v
       if i > 1
           # Solve H[1:i+1,1:i]*x = bnorm*e_1, using least squares
           # TODO: Implement Givens rotations
           rhs = zeros(Float64,i+1)
           rhs[1] = bnorm
           x = H[1:i+1,1:i]\rhs
           res = norm(H[1:i+1,1:i]*x-rhs)
           print("iteration = ")
           print(i)
           print(", residual = ")
           println(res)
           if res < tol
               return [Q,x]
           end
       end
    end
    println("GMRES did not terminate")
    return[Q,x]
end

function GMRES_verbose(A,b,inner,tol,n,cond)

    nom = a -> sqrt(abs(inner(a,a)))
    H = zeros(Complex{Float64},n+1,n)
    bnorm = nom(b)
    x = 0.
    Q = [(1.0/bnorm)*b]
    for i = 1:n
       #tic()
       println("Operator application: ")
       @time v = A(Q[i])
       #toc()
       #tic()
       println("Inner products: ")
       for j = 1:i
           #tic()
           @time H[j,i] = inner(Q[j],v)
           #toc()
           v = cond(v - H[j,i]*Q[j])
       end
       v = cond(v)
       println("Assembling Q:")
       @time H[i+1,i] = nom(v)
       @time Q = vcat(Q,[copy((1.0/H[i+1,i])*v)])
       #print("Arnoldi: ")
       #toc()
       #return v
       if i > 1
           # Solve H[1:i+1,1:i]*x = bnorm*e_1, using least squares
           # TODO: Implement Givens rotations
           rhs = zeros(Float64,i+1)
           rhs[1] = bnorm
           println("Solve least squares")
           @time x = H[1:i+1,1:i]\rhs
           res = norm(H[1:i+1,1:i]*x-rhs)
           print("iteration = ")
           print(i)
           print(", residual = ")
           println(res)
           if res < tol
               return [Q,x]
           end
       end
    end
    println("GMRES did not terminate")
    return[Q,x]
end

end # module
