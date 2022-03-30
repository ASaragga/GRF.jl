# Modelos de Volatilidade
""" 
    EWMA(retornos, lambda)  

Calcula médias móveis exponencialmente ponderadas

## Argumentos
* `retornos`: retornos da carteira
* `lambda`: parâmetro de alisamento exponencial
"""
function EWMA(retornos, lambda)
    m = length(retornos)
    media = mean(retornos)                 
    sigma2 = 0
    for k in 1:m
        sigma2 += (1-lambda)*lambda^(k-1) * (retornos[end-k+1] - media)^2
    end
    return sqrt(sigma2)                     
end


# Abordagem Empírica
""" 
    VaR(retornos, alfa, V)  

Calcula o Valor-em-Risco (VaR)

## Argumentos
* `retornos`: retornos da carteira
* `alfa`: nível de significância
* `V`: valor da carteira
"""
VaR(retornos,alfa,V) = -quantile(retornos,alfa) * V


""" 
    ETL(retornos, alfa, V)  

Calcula a Perda Esperada na Cauda (ETL)

## Argumentos
* `retornos`: retornos da carteira
* `alfa`: nível de significância
* `V`: valor da carteira
"""
function ETL(retornos, alfa, V)
    nobs = length(retornos)
    corte = quantile(retornos, alfa)
    excedimentos = 0
    n_excedimentos = 0
    for i in 1:nobs
        if retornos[i] < corte
            excedimentos += retornos[i]
            n_excedimentos += 1
        end
    end
    return - (excedimentos/n_excedimentos) * V
end


# Método Analítico
""" 
    VaR0(alfa, mu, sigma, V)

Calcular o VaR-zero

## Argumentos
* `alfa`: nível de significância
* `mu`: retorno esperado da carteira
* `sigma`: desvio-padrão dos retornos da carteira
* `V`: valor da carteira
"""
VaR0(alfa, mu, sigma,V) = (-quantile(Normal(0,1), alfa) * sigma - mu) * V


""" 
    VaRm(alfa, sigma, V)

Calcular o VaR-média

## Argumentos
* `alfa`: nível de significância
* `sigma`: desvio-padrão dos retornos da carteira
* `V`: valor da carteira
"""
VaRm(alfa, sigma, V) = -quantile(Normal(0,1), alfa) * sigma * V


""" 
    ETL(alfa, mu, sigma, V)

Calcular a Perda Esperda na Cauda

## Argumentos
* `alfa`: nível de significância
* `mu`: retorno esperado da carteira
* `sigma`: desvio-padrão dos retornos da carteira
* `V`: valor da carteira
"""
ETL(alfa, mu, sigma, V) = (mu + sigma * pdf(Normal(0,1), quantile(Normal(0,1), alfa)) / alfa) * V


""" 
    semaforo(retornos, alfa)

Teste do semáforo de Basileia
""" 
function semaforo(retornos, alfa)
    nobs = length(retornos)
    excep = zeros(nobs)
    sigma = std(retornos)
    RaR = - sigma * quantile(Normal(0,1), alfa)
    for i in 1:nobs
        if retornos[i] <= -RaR
            excep[i] = 1
        end
    end
    return sum(excep)
end
