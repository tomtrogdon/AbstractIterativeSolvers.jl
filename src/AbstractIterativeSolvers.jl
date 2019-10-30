module AbstractIterativeSolvers

using LinearAlgebra

export GMRES

function GMRES(A,b,inner,tol,n)
    nom = a -> sqrt(abs(inner(a,a)))
    H = zeros(Complex{Float64},n+1,n)
    bnorm = nom(b)
    x = 0.
    Q = [(1.0/bnorm)*b]
    for i = 1:n
       #tic()
       v = A(Q[i])
       #print("Operator application: ")
       #toc()
       #tic()
       for j = 1:i
           #tic()
           H[j,i] = inner(Q[j],v)
           #toc()
           v = v - H[j,i]*Q[j]
       end
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

end # module
