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


function cvm(Sigma)
    n = size(Sigma,1)
    modelo = Model(OSQP.Optimizer)
    set_silent(modelo)  
    @variable(modelo, w[1:n])                      
    @objective(modelo, Min, w' * Sigma * w) 
    @constraint(modelo, sum(w) == 1)        
    @constraint(modelo, w .>= 0)            
    optimize!(modelo)                       
    return value.(w)
end


function fe(mu, Sigma, muk)
    n = size(Sigma,1)
    modelo = Model(OSQP.Optimizer)
    set_silent(modelo)  
    @variable(modelo, w[1:n])     
    @objective(modelo, Min, w' * Sigma * w)  
    @constraint(modelo, sum(w) == 1)  
    @constraint(modelo, w .>= 0)      
    @constraint(modelo, mu' * w >= muk)
    optimize!(modelo)                   
    return value.(w)                                 
  end



function gfe(mu, Sigma, ncarteiras = 10) 
    muk = zeros(ncarteiras) 
    vark = zeros(ncarteiras)                        
    mu_max = maximum(mu) 

    w_cvm = cvm(Sigma) 
    muk[1] = mu' * w_cvm  # MVP: mu
    vark[1] = w_cvm' * Sigma * w_cvm  # MVP: sigma^2
    inc = (mu_max-muk[1])/(ncarteiras-1)                              
    for i in 2:ncarteiras
        muk[i] = muk[1] + inc * (i-1)
        wk = fe(mu, Sigma, muk[i])                               
        vark[i] = wk' * Sigma * wk 
    end
    vars = diag(Sigma)
    fig = plot(vark, muk, xlabel = L"Risco ($\sigma^2$)", ylabel = "Valor Esperado do Retorno", label = "Fronteira Eficiente", xlim = (0, maximum(vars) * 1.1), ylim = (0, maximum(mu)*1.1), framestyle = :box, legend = :bottomright)
    fig = scatter!(vars, mu, label = "Ativos")
    fig = scatter!([vark[1]], [muk[1]], label = "Carteira variância minima")
    return fig
end

function alocar(mu, Sigma, lista, ncarteiras = 10)
    modelo = Model(OSQP.Optimizer)
    set_silent(modelo)       
    n = length(mu)                           
    @variable(modelo, w[1:n])                           
    @objective(modelo, Min, w' * Sigma * w)             
    @constraint(modelo, sum(w) == 1)                    
    @constraint(modelo, w .>= 0)                        
    optimize!(modelo)                                   
    w_CVM = value.(w)                                   
    mu_k = zeros(ncarteiras)
    w_k = zeros(n,ncarteiras)
    mu_CVM = mu' * w_CVM
    for i in 0:ncarteiras-1
        modelo = Model(OSQP.Optimizer)
        set_silent(modelo)                                  
        @variable(modelo, w[1:n])                           
        @objective(modelo, Min, w' * Sigma * w)             
        @constraint(modelo, sum(w) == 1)                    
        @constraint(modelo, w .>= 0)                        
        mu_max = maximum(mu)
        incremento = (mu_max-mu_CVM)/(ncarteiras-1)
        mu_k[i+1] = mu_CVM + incremento * i
        @constraint(modelo, mu' * w >= mu_k[i+1])                
        optimize!(modelo)                                  
        w_k[:,i+1] = value.(w)   
    end
    ticklabel = "P" .* string.(collect(1:ncarteiras))
    ticklabel[1] = "CVM"
    ticklabel[ncarteiras] = "CRM"
    fig = groupedbar(w_k', bar_position = :stack, bar_width=1.0, xlabel = "Carteiras", xticks=(1:ncarteiras, ticklabel), ylabel = "Pesos", label = permutedims(String.(lista)), legend = :bottomright, framestyle = :box)
    return fig
end


#=
function installedx()
    @warn "Pkg.installed() is deprecated"
    deps = dependencies()
    installs = Dict{String, VersionNumber}()
    for (uuid, dep) in deps
        dep.is_direct_dep || continue
        dep.version === nothing && continue
        installs[dep.name] = dep.version
    end
    return installs
end
=#
