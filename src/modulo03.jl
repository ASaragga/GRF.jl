function ETL(retornos, alfa, V)
    nobs = length(retornos)
    corte = quantile(retornos, alfa)
    excedimentos = 0
    n_excedimentos = 0
    for i in 1:nobs
        if retornos[i] <= corte
            excedimentos += retornos[i]
            n_excedimentos += 1
        end
    end
    return - (excedimentos/n_excedimentos - mean(retornos)) * V
end


function EWMA(retornos, lambda)
    m = length(retornos)
    media = mean(retornos)                 
    sigma2 = 0
    for k in 1:m
        sigma2 += (1-lambda)*lambda^(k-1) * (retornos[end-k+1] - media)^2
    end
    return sqrt(sigma2)                     
end


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
