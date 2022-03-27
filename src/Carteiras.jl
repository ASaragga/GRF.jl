# Dados
using CSV
using DataFrames       
# Gráficos
using StatsPlots       
# Estatistica
using GLM
using Statistics       
# Optimização e Álgebra Linear
using Ipopt            
using HiGHS
using JuMP             
using LinearAlgebra   

using GRF

DJ30r = CSV.read("/Users/antonio/Documents/Universidade/GRF-2022/Dados/DJ30r.csv", DataFrame);
carteira = Matrix(DJ30r[!,Not([:SP500, :Datas])]);

beta = zeros(30);
beta = [coef(lm(term(i) ~ term(:SP500), DJ30r))[2] for i in names(DJ30r[!, Not([:SP500, :Datas])])];  # Vetor de Betas
rf = 0
rm = 0.05
mu = rf .+ beta * (rm - rf);  # Vetor de retornos esperados de acordo com CAPM

Sigma = cov(carteira);    # Matriz de covariâncias 
w_cvm = GRF.cef(Sigma)    # Carteira de variância mínima
r_cvm =  w_cvm' * mu      # Valor esperado do retorno da carteira de variância mínima
sigma_cvm = sqrt( w_cvm' * Sigma * w_cvm)  # DP retorno da carteira de variância mínima
println("Volatilidade(CVM) = ", sigma_cvm, ", Retorno esperado(CVM) = ", r_cvm)

fig1 = GRF.gfe(Sigma, mu)  # Fronteira eficiente
display(fig1)
fig2 = GRF.alocar(Sigma, mu, names(DJ30r[!, Not([:SP500, :Datas])]))
display(fig2)
