
using Plots; pgfplots()
using Random
using LaTeXStrings
using LinearAlgebra
using NLsolve              # To solve non linear systems
using Distributions        # For markov process function
using SpecialFunctions     # For Tauchen algorithm function
using Dierckx              # Spline interpolation


# PARAMETERS
# ==========

global chi  = 1.5          # Scale parameter of labor
global psi  = 1.0          # Elasticity of labor
global rho  = 1.01         # Elasiticty of consumption
global beta = 0.9

phi_grid    = range(0,     # Lagrange multiplier of the implementability constraint
              stop=1,
              length=30)     

global G    = [0.1 0.2]'   # State-space of government expenditures

global Pi   = [0.9 0.1;    # Markov transition matrix
               0.5 0.5];

# ======================
# 1.a Defining functions
# ======================

function tau(phi)
    1 .- (1 .+ phi .* (1 - rho)) ./ (1 .+ phi .* (1 + psi))
end

# =======================
# 1.a Solving the problem
# =======================

plot(phi_grid, tau(phi_grid),
     title="Tax rate as function of phi",
     label="",
     xlabel=L"$\phi$",
     ylabel=L"$\tau(\phi)$")

# ======================
# 1.b Defining functions
# ======================
"""
    starEquation(c, n, phi, b0)

Return the value of the star equation for a given consumption `c`, labor `n`, lagrange multiplier `phi`,
and initial debt b0. For t>0, use b0=0.

"""
function starEquation(c, n, phi, b0::Any=0)
    Uc  = c ^ (- rho)
    Ucc = - rho * c ^ (-rho - 1)
    Un  = - chi * n ^ psi
    Unn = - chi * psi * n ^ (psi - 1)
    
    return (1 + phi) * (Uc + Un) + phi * Ucc * (c - b0) + phi * Unn * n
end

"""
    timeInvariantAllocation(g, phi)

Solves the time invariant allocation for a given government expenditure g and lagrange multiplier phi.
"""
function timeInvariantAllocation(g, phi)

    function f!(F, x)
        F[1] = x[1] - x[2] + g
        F[2] = starEquation(x[1], x[2], phi, 0)
    end

    nlsolve(f!, [0.5; 0.5])
end

"""
Compute the time invariant allocation for a given set of governement expenditures and lagrange multipliers

Arguments:
    G: a vector of possible values for the government expenditures
    phi: a vector of possible values for phi

Return:
    (C, N): a tuple of matrices where X[i, j] is the time invariant allocation of X for g[i] and phi[j]
"""
function timeInvariantAllocationVector(G, phi)
    C = zeros(length(G), length(phi))
    N = zeros(length(G), length(phi))

    for (i, g) in enumerate(G)
        for (j, phi) in enumerate(phi)
            roots   = timeInvariantAllocation(g, phi).zero
            C[i, j] = roots[1]
            N[i, j] = roots[2]
        end
    end

    return (C', N')
end

# =======================
# 1.b Solving the problem
# =======================
C, N = timeInvariantAllocationVector(G, phi_grid)

plot(
    layout=(2,1),
    plot(phi_grid, C,
        label=[L"g_L" L"g_H"],
        ylabel=L"$c(\Phi)$",
        title="Time invariant allocation"
    ),
    plot(phi_grid, N,
        label=[L"g_L" L"g_H"],
        ylabel=L"$n(\Phi)$",
        xlabel=L"$\Phi$"
    )
)

# ======================
# 1.c Defining Functions
# ======================
"""
    marginalUtilities(Consumption, Labor)
    
"""
function marginalUtilities(c::Adjoint, n::Adjoint)
    Uc = c.^(-rho)
    Un = -chi .* n .^ psi
    
    return (Uc, Un)
end

"""
    bonds(G, phi, Uc, Un, Pi, beta)

Returns:
    B : a (Nphi x NStates) matrix where B(i,j) is b(phi[i], G[i])
"""
function bonds(G, Pi, phi, beta)
    Nstates = length(G)
    Nphi    = length(phi)
    C, N    = timeInvariantAllocationVector(G, phi)
    Uc, Un  = marginalUtilities(C,N)
    B  = zeros(Nstates, Nphi)
    Id = Matrix{Float64}(I, Nstates, Nstates)
    for i in 1:Nphi #Iterate over all values of phi
        c = C[i,:]
        n = N[i,:]
        uc = Uc[i,:]
        un = Un[i,:]
        
        # The system writes:
        # b = Ab + K
        A  = beta .* (1 ./ uc) * uc' .* Pi
        k  = c .+ un./uc .* n
        
        B[:,i] = inv(Id - A) * k # We solve the system for a given phi
    end
    
    return B'
end

# =======================
# 1.c Solving the problem
# =======================
B = bonds(G, Pi, phi_grid, beta);

plot(phi_grid, B,
    xlabel=L"$\Phi$",
    ylabel=L"$b(\Phi)$",
    label=[L"$g_L$", L"$g_H$"],
    title="Bonds in function of lagrange multiplier")

# ======================
# 2.a Defining functions
# ======================
"""
Solves the time invariant allocation for a given government expenditure g and lagrange multiplier phi
"""
function timeZeroAllocation(g, phi, b0)

    function f!(F, x)
        F[1] = x[1] - x[2] + g 
        F[2] = starEquation(x[1], x[2], phi, b0)
    end

    nlsolve(f!, [0.5; 0.5])
end

"""
Solves the time zero allocation for a vector of phi
"""
function timeZeroAllocationVector(g, phi, b0)
    C0 = zeros(length(phi), length(b0))
    N0 = zeros(length(phi), length(b0))
    
    for (i, phi) in enumerate(phi)
        for (j, b0) in enumerate(b0)
            roots   = timeZeroAllocation(g, phi, b0).zero
            C0[i, j] = roots[1]
            N0[i, j] = roots[2]
        end
    end
    
    return (C0, N0)
end

# =======================
# 2.a Solving the problem
# =======================
b0 = range(-0.1, stop=0.1, length=length(phi_grid))
C0, N0 = timeZeroAllocationVector(G[1], phi_grid, b0)

levels = collect(0.35:0.04:0.9)
plot(layout=(1,2),
    plot(b0, phi_grid, C0',
        linetype=:contour,
        xlabel= L"b_0",
        ylabel= L"\phi",
        title=L"$C_{0_{\phi, b_0}}$",
        levels = levels,
        legend=false
    ),
    plot(b0, phi_grid, N0',
        linetype=:contour,
        xlabel= L"b_0",
        title=L"$N_{0_{\phi, b_0}}$",
        levels = levels,
    )
)

# ======================
# 2.b Defining functions
# ======================
"""
Bisection algorithm for root finding
"""
function bisec(f::Function, a, b, tol::AbstractFloat=1e-5, maxiter::Integer=100)
    fa = f(a)
    fa*f(b) <= 0 || error("No real root in [a,b]")
    i = 0
    local c
    while b-a > tol
        i += 1
        i != maxiter || error("Max iteration exceeded")
        c = (a+b)/2
        fc = f(c)
        if fc == 0
            break
        elseif fa*fc > 0
            a = c  # Root is in the right half of [a,b].
            fa = fc
        else
            b = c  # Root is in the left half of [a,b].
        end
    end
    return c
end

"""
Compute the implementability constraint for the optimal allocation (c*, n*)
"""
function implementabilityConstraint(phi, b0)
    c, n = timeZeroAllocation(G[1], phi, b0).zero
    Uc  = c ^ (- rho)
    Ucc = - rho * c ^ (-rho - 1)
    Un  = - chi * n ^ psi
    Unn = - chi * psi * n ^ (psi - 1)
    return Uc * c + Un * n - Uc * b0
end

"""
Returns the equilibrium value of phi given an initial level of debt b0
"""
function phiFromb0(b0)
    phi_true = bisec(x -> implementabilityConstraint(x, b0), 0, 1)
    return phi_true
end

# =======================
# 2.b Solving the problem
# =======================
# define a grid for b_0 on ]-0.1, 0.1] (there is no solution when b_0 = -0.1)
b0_grid = range(-0.099, stop = 0.1, length = 20)
phi_true_grid = zeros(length(b0_grid)) # phi(b_0) vector preallocation

# Solving the Implementability constraint for each values of b_0 using a bisection algorithm
for (i, b0) in enumerate(b0_grid)
    phi_true_grid[i] = phiFromb0(b0)
end

plot(b0_grid, phi_true_grid,
    title="Lagrange Multiplier as Function of Initial Debt",
    xlabel=L"$b_0$",
    ylabel=L"\phi(b_0)",
    label="")

# ====================
# 3 Defining functions
# ====================
function markovChainSimulation(transitionMatrix::Array, length::Integer, seed::Integer=42, initialState::Integer=1)
    chain = Array{Int32,1}(undef, length)
    chain[1] = initialState
    Nstates  = size(transitionMatrix)[1]
    Random.seed!(seed)
    for i in 2:length
        distribution = Multinomial(1, transitionMatrix[chain[i-1],:])
        chain[i]     = rand(distribution)' * (1:Nstates)
    end

    return chain
end

function simulateEconomy(states, b0::Float64, transitionMatrix::Array, N::Integer, seed::Integer=42, initialState::Integer=1)
    g, c, b, t, n = [zeros(N) for i in 1:6]
    
    chain = markovChainSimulation(transitionMatrix, N, seed, initialState)
    chain[1] = 0 # To avoid time 0 overwriting
    
    # Allocations in t=0
    # ==================
    g[1] = states[initialState]
    b[1] = b0
    phi  = phiFromb0(b0)
    c[1], n[1] = timeZeroAllocation(g[1], phi, b[1]).zero
    t[1] = (c[1]^(-rho) - chi * n[1]^psi) / (c[1]^(-rho))
    t[2:end] .= tau(phi) # Tax rate is constant in t>1
    
    # Compute allocations for each states
    # ===================================
    for (i,) in enumerate(states)
        g[chain .== i] .= states[i]
        cons, labor = timeInvariantAllocation(states[i], phi).zero
        c[chain .== i] .= cons
        n[chain .== i] .= labor
        
        # We use spline approx to compute bonds from our previous computed matrix
        splineBonds = Spline1D(phi_grid, B[:,i])
        b[chain .== i] .= splineBonds(phi)
    end
    
    return g, c, n, b, t
end

function plotSimulatedEconomy(states, b0, transitionMatrix, periods, seed::Integer=42)
    g, c, n, b, t = simulateEconomy(states, b0, transitionMatrix, periods, seed);

    plot(layout=(3,2),
        plot(c, title="Consumption", label=""),
        plot(n, title="Labor supply", label=""),
        plot(b, title="Government debt", label=""),
        plot(round.(t, digits=2), title="Tax Rate", label=""),
        plot(g, title="Government Spending", label=""),
        plot(b./n, title="Debt/output Ratio", label=""),
    )
end;

plotSimulatedEconomy(G, 0.0, Pi, 100)

plotSimulatedEconomy(G, 0.1, Pi, 100)

plotSimulatedEconomy(G, -0.09999, Pi, 100)

# ======================
# 3.d Defining Functions
# ======================

"""
Approximate AR1 with finite markov process
Proudly stolen on quantEcon github repository
https://github.com/QuantEcon/QuantEcon.jl/blob/master/src/markov/markov_approx.jl
"""

std_norm_cdf(x::T) where {T <: Real} = 0.5 * erfc(-x/sqrt(2))
std_norm_cdf(x::Array{T}) where {T <: Real} = 0.5 .* erfc(-x./sqrt(2))

"""
Tauchen's (1996) method for approximating AR(1) process with finite markov chain
The process follows
    y_t = mu + rho y_{t-1} + epsilon_t

where epsilon_t sim N (0, sigma^2)

    Args:
        - N::Integer: Number of points in markov process
        - ρ::Real : Persistence parameter in AR(1) process
        - σ::Real : Standard deviation of random component of AR(1) process
        - μ::Real(0.0) : Mean of AR(1) process
        - n_std::Integer(3) : The number of standard deviations to each side the process
          should span
"""
function tauchen(N::Integer, ρ::Real, σ::Real, μ::Real=0.0, n_std::Integer=3)
    # Get discretized space
    a_bar = n_std * sqrt(σ^2 / (1 - ρ^2))
    y = range(-a_bar, stop=a_bar, length=N)
    d = y[2] - y[1]

    # Get transition probabilities
    Π = zeros(N, N)
    for row = 1:N
        # Do end points first
        Π[row, 1] = std_norm_cdf((y[1] - ρ*y[row] + d/2) / σ)
        Π[row, N] = 1 - std_norm_cdf((y[N] - ρ*y[row] - d/2) / σ)

        # fill in the middle columns
        for col = 2:N-1
            Π[row, col] = (std_norm_cdf((y[col] - ρ*y[row] + d/2) / σ) -
                           std_norm_cdf((y[col] - ρ*y[row] - d/2) / σ))
        end
    end

    # NOTE: I need to shift this vector after finding probabilities
    #       because when finding the probabilities I use a function
    #       std_norm_cdf that assumes its input argument is distributed
    #       N(0, 1). After adding the mean E[y] is no longer 0, so
    #       I would be passing elements with the wrong distribution.
    #
    #       It is ok to do after the fact because adding this constant to each
    #       term effectively shifts the entire distribution. Because the
    #       normal distribution is symmetric and we just care about relative
    #       distances between points, the probabilities will be the same.
    #
    #       I could have shifted it before, but then I would need to evaluate
    #       the cdf with a function that allows the distribution of input
    #       arguments to be [μ/(1 - ρ), 1] instead of [0, 1]

    yy = y .+ μ / (1 - ρ) # center process around its mean (wbar / (1 - rho)) in new variable

    # renormalize. In some test cases the rows sum to something that is 2e-15
    # away from 1.0, which caused problems in the MarkovChain constructor
    Π = Π./sum(Π, dims = 2)

    return Π, yy
end

# =======================
# 3.d Solving the problem
# =======================
Pi_large, G_large = tauchen(100, 0.95, 0.08, 0.1)
G_large = G_large./10
B = bonds(G_large, Pi_large, phi_grid, beta)

plotSimulatedEconomy(G_large, 0.1, Pi_large, 100)
