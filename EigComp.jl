using LinearAlgebra

function eigcomp(A)
    n = size(A, 1)
    u = ones(n, 1)
    y = ones(n, 1)
    beta1 = 0

    eta = norm(u)
    y = u / eta
    u = A * y
    beta2 = dot(y, u)

    while abs((beta2 - beta1) / beta1) > 1e-12
        beta1 = beta2
        eta = norm(u)
        y = u / eta
        u = A * y
        beta2 = dot(y, u)
    end

    return beta2, y
end
