using Plots, Random, HypothesisTests, Statistics, LinearAlgebra, Distributions, StatsBase, Optim, Plots, LaTeXStrings, KernelDensity, DataFrames, GLM, StatsPlots
Random.seed!(6413)
################################
#                              #
#          Problem 1           #
#                              #
################################

println("PROBLEM 1")
# generate random numbers from a generalized lambda distribution given a uniform distribution
gld(y, lambda1, lambda2, lambda3, lambda4) = lambda1 + (y^lambda3 - (1 - y)^lambda4) / lambda2

# test statistic function for the Jarque Bera test
function test_statisc(x)
    length(x) / 6 * (skewness(x)^2 + 1 / 4 * kurtosis(x)^2)
end

##
B = 1000 # number of monte carlo simulations
T = 200 # length of each sequence
lambda1 = 0
lambda2 = 0.197454
lambda3 = 0.134915
lambda4 = 0.134915
level = 0.05


# sample from uniform distribution
rand_uniform = rand(B, T)

gld_samples = gld.(rand_uniform, lambda1, lambda2, lambda3, lambda4)
stats = test_statisc.(eachrow(gld_samples))

test_jb = ccdf.(Chisq(2), stats) .<= level # reject if probability is lower than threshold
println("Size for T=200: ", mean(test_jb) * 100)

## repeat for T=500
T = 500
rand_uniform = rand(B, T)
gld_samples = gld.(rand_uniform, lambda1, lambda2, lambda3, lambda4)
stats = test_statisc.(eachrow(gld_samples))
test_jb = ccdf.(Chisq(2), stats) .<= level
println("Size for T=500: ", mean(test_jb) * 100)

## repeats for assymetric fat tailed distribution

B = 1000 # number of monte carlo simulations
T = 200 # length of each sequence
lambda1 = 0
lambda2 = 1
lambda3 = 0.000070
lambda4 = 0.10000
level = 0.05
rand_uniform = rand(B, T)

gld_samples = gld.(rand_uniform, lambda1, lambda2, lambda3, lambda4)
stats = test_statisc.(eachrow(gld_samples))

test_jb = ccdf.(Chisq(2), stats) .<= level # reject if probability is lower than threshold
mean(test_jb) * 100
println("Power for T=200: ", mean(test_jb) * 100)

## repeat for T=500
T = 500 # length of each sequence
lambda1 = 0
lambda2 = 1
lambda3 = 0.000070
lambda4 = 0.10000
level = 0.05
rand_uniform = rand(B, T)
gld_samples = gld.(rand_uniform, lambda1, lambda2, lambda3, lambda4)
stats = test_statisc.(eachrow(gld_samples))
test_jb = ccdf.(Chisq(2), stats) .<= level # reject if probability is lower than threshold
println("Power for T=500: ", mean(test_jb) * 100)
##

################################
#                              #
#          Problem 2           #
#                              #
################################
println("")
println("PROBLEM 2")

function gbar(θ, x)
    μ, σ_2 = θ
    return [mean(x .- μ), mean((x .- μ) .^ 2 .- σ_2), mean((x .- μ) .^ 3), mean((x .- μ) .^ 4 .- 3 * σ_2^2)]
end

function obj_function_gmm(θ, x)
    μ, σ_2 = θ
    S = [σ_2 0 3σ_2^2 0; 0 2σ_2^2 0 12σ_2^3; 3σ_2^2 0 15σ_2^3 0; 0 12σ_2^3 0 96σ_2^4]
    return gbar(θ, x)' * inv(S) * gbar(θ, x) * length(x)
end

function gradient(f, point, eps1, eps2)
    x, y = point
    g1 = (f([x + eps1, y]) - f([x - eps1, y])) / (2 * eps1)
    g2 = (f([x, y + eps2]) - f([x, y - eps2])) / (2 * eps2)
    return [g1, g2]
end

function hessian(f, point, eps1, eps2)
    x, y = point
    h11 = (f([x + eps1, y]) - 2 * f([x, y]) + f([x - eps1, y])) / (eps1^2)
    h12 = (f([x + eps1, y + eps2]) - f([x + eps1, y - eps2]) - f([x - eps1, y + eps2]) + f([x - eps1, y - eps2])) / (4 * eps1 * eps2)
    h22 = (f([x, y + eps2]) - 2 * f([x, y]) + f([x, y - eps2])) / (eps2^2)
    return [h11 h12; h12 h22]
    return

end

function newton_manual(f, θ0, stop_criter)
    eps1 = 1e-5
    eps2 = 1e-5
    iter_max = 1000
    iter = 1
    err = stop_criter * 2
    while iter <= iter_max && err > stop_criter
        g = gradient(f, θ0, eps1, eps2)
        h = hessian(f, θ0, eps1, eps2)
        θ_new = θ0 .- inv(h) * g

        err = maximum(abs.(θ_new - θ0))
        θ0 = θ_new
        iter = iter + 1
    end
    return θ0
end

N = 10000
mu = 2
sigma = 4
x = rand(Normal(mu, sigma), N)
θ0_optim = optimize(θ -> obj_function_gmm(θ, x), [1.1, 1.1]).minimizer
θ0_newton = newton_manual(θ -> obj_function_gmm(θ, x), [0.0, 4], 1e-6)
println("Results for Optim ")
println(θ0_optim)

println("Results for NR own algorithm:")
println(θ0_newton)

##
################################
#                              #
#          Problem 3           #
#                              #
################################

println("")
println("PROBLEM 3")

alpha = 0.5
beta = 0.8
sigma = 1
N = 1000

function simulate_data(alpha, beta, sigma, N)
    inov = rand(Normal(0, sigma), N)
    y = zeros(N)
    inov[1] = 0
    for i in 2:length(y)
        y[i] = alpha * y[i-1] + inov[i] + beta * inov[i-1]
    end
    return y
end


function theoretical_covs(α, β, σ)
    gamma = zeros(3)
    gamma[1] = σ^2 * (2 * α * β + 1 + β^2) / (1 - α^2)
    gamma[2] = α * σ^2 * (2 * α * β + 1 + β^2) / (1 - α^2) + β * σ^2
    gamma[3] = α^2 * σ^2 * (2 * α * β + 1 + β^2) / (1 - α^2) + α * β * σ^2
    return gamma

end

function cov_estimation_obj(α, β, σ, gamma_hat, W)
    m = zeros(3)
    gammas = theoretical_covs(α, β, σ)
    m[1] = gamma_hat[1] - gammas[1]
    m[2] = gamma_hat[2] - gammas[2]
    m[3] = gamma_hat[3] - gammas[3]
    return m' * W * m
end

function resids(alpha, beta, y, e_0)
    e = zeros(length(y))
    e[1] = e_0
    for i in 2:length(y)
        e[i] = y[i] - alpha * y[i-1] - beta * e[i-1]
    end
    return e
end

## part i
y = simulate_data(alpha, beta, sigma, N)

function obj_function_nls(param)
    res = resids(param[1], param[2], y, 0)
    return res' * res
end
println("NLS estimation")
sol = optimize(obj_function_nls, [0.2, 0.2]).minimizer
# given alpha and beta find the std or res
sigma2 = var(resids(sol[1], sol[2], y, 0))
println("alpha = ", sol[1])
println("beta = ", sol[2])
println("sigma2 = ", sigma2)


## part ii

function gmm_part2(alpha, beta, sigma, y, W)
    e = resids(alpha, beta, y, 0)
    moments = [var(e) - sigma^2, autocov(e, [1]), autocov(e, [2])]
    return moments' * W * moments
end
println("")
println("GMM estimation")
sol = optimize(param -> gmm_part2(param[1], param[2], param[3], y, I), [0.4, 0.0, 3]).minimizer
println("alpha = ", sol[1])
println("beta = ", sol[2])
println("sigma2 = ", sol[3])
## part iii


B = 1000
sols = zeros(B, 3)
for j in 1:B
    inov = rand(Normal(0, sigma), N)
    local y = simulate_data(alpha, beta, sigma, N)
    gammas = autocov(y, [0, 1, 2])
    cov_aux(param) = cov_estimation_obj(param[1], param[2], param[3], gammas, I)
    sols[j, :] = optimize(cov_aux, [0.5, 0, 4]).minimizer
end
println("")
println("Covariance estimation")
println("alpha = ", sols[1, 1])
println("beta = ", sols[1, 2])
println("sigma2 = ", sols[1, 3])

alphas = sols[:, 1]
histogram(alphas, normalize=true, label="Data")
plot!(fit(Normal, alphas), label="Normal distribution")
xlabel!(L"\hat \alpha")
savefig("alpha_hat_density.pdf")
##

beta = -0.45
B = 1000
sols = zeros(B, 3)
for j in 1:B
    inov = rand(Normal(0, sigma), N)
    local y = simulate_data(alpha, beta, sigma, N)
    gammas = autocov(y, [0, 1, 2])
    cov_aux(param) = cov_estimation_obj(param[1], param[2], param[3], gammas, I)
    sols[j, :] = optimize(cov_aux, [0.5, 0, 4]).minimizer
end
alphas = sols[:, 1]
histogram(alphas, normalize=true, label="Data")
plot!(fit(Normal, alphas), label="Normal distribution")
xlabel!(L"\hat \alpha")
savefig("alpha_hat_density_beta_5.pdf")

alpha_mean = mean(sols[:, 1])
beta_mean = mean(sols[:, 2])
sigma_mean = mean(sols[:, 3])
theoretical_covs(alpha_mean, beta_mean, sigma_mean)

################################
#                              #
#          PROBLEM 4           #
#                              #
################################

println("")
println("PROBLEM 4")

T = 1000
σ_e = 1
θ = 0.5
inov = rand(Normal(0, σ_e), T)
y = zeros(T)
inov[1] = 0
for i in 2:length(y)
    y[i] = inov[i] + θ * inov[i-1]
end

function lag(vector, n)
    return [ones(n) * NaN; vector[1:end-n]]

end

function obj_function_md(param, ϕ_1, ϕ_2, σ_u, W)
    θ, σ_e = param
    phi_1_theory = θ * (θ^2 + 1) / (θ^4 + θ^2 + 1)
    phi_2_theory = -θ^2 / (θ^4 + θ^2 + 1)
    σ2_u_theory = σ_e^2 * (θ^6 + θ^4 + θ^2 + 1) / (θ^4 + θ^2 + 1)
    moment = [ϕ_1 - phi_1_theory, ϕ_2 - phi_2_theory, σ_u^2 - σ2_u_theory]
    return moment' * W * moment
end


y_1 = lag(y, 1)
y_2 = lag(y, 2)
data = DataFrame(y=y[3:end], y_1=y_1[3:end], y_2=y_2[3:end])
ols = lm(@formula(y ~ -1 + y_1 + y_2), data)
println("OLS estimation")
println(ols)
var_u = var(residuals(ols))
println("Sigma2 u= ", var_u)

ϕ_1 = coef(ols)[1]
ϕ_2 = coef(ols)[2]
σ_u = sqrt(var_u)
println("")
sol = optimize(param -> obj_function_md(param, ϕ_1, ϕ_2, σ_u, I), [0.8, 2]).minimizer
println("Minimum distance estimation:")
println("theta = ", sol[1])
println("sigma_e = ", sol[2])
