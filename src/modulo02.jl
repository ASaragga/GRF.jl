function cef(Sigma, mu, muk)  # Carteira efficiente para (muk)
    n = size(Sigma,1)
    modelo = Model(Ipopt.Optimizer)
    set_silent(modelo)  
    @variable(modelo, w[1:n])     
    @objective(modelo, Min, w' * Sigma * w)  
    @constraint(modelo, sum(w) == 1)  
    @constraint(modelo, w .>= 0)      
    @constraint(modelo, mu' * w >= muk)
    optimize!(modelo)                   
    return value.(w)                                 
end

function cef(Sigma)   # Carteira variância mínima
    n = size(Sigma,1)
    modelo = Model(Ipopt.Optimizer)
    set_silent(modelo)  
    @variable(modelo, w[1:n])                      
    @objective(modelo, Min, w' * Sigma * w) 
    @constraint(modelo, sum(w) == 1)        
    @constraint(modelo, w .>= 0)            
    optimize!(modelo)                       
    return value.(w)
end

function fef(Sigma, mu, ncarteiras = 10) # Fronteira eficiente
    muk = zeros(ncarteiras) 
    vark = zeros(ncarteiras)
    wk = zeros(size(Sigma,1), ncarteiras)
    mu_max = maximum(mu) 
    wk[:,1] = cef(Sigma) 
    muk[1] = mu' * wk[:,1]  # CVM: mu
    vark[1] = wk[:,1]' * Sigma * wk[:,1]  # CVM: sigma^2
    inc = (mu_max-muk[1])/(ncarteiras-1)                              
    for i in 2:ncarteiras
        muk[i] = muk[1] + inc * (i-1)
        wk[:,i] = cef(Sigma, mu, muk[i])                               
        vark[i] = wk[:,i]' * Sigma * wk[:,i]
    end
    return wk, vark, muk
end

function gfe(Sigma, mu, ncarteiras = 10) 
    dpfe = sqrt.(252 * fef(Sigma, mu, ncarteiras)[2])
    mufe = fef(Sigma, mu, ncarteiras)[3]
    DPs = sqrt.(252 * diag(Sigma))
    fig = plot(dpfe, mufe, xlabel = "DP[r] Anualizado", ylabel = L"\mathbb{E}[r]", label = "Fronteira Eficiente", xlim = (0, maximum(DPs) * 1.1), ylim = (0, maximum(mu)*1.1), framestyle = :box, legend = :bottomright)
    fig = scatter!(DPs, mu, label = "Ativos")
    fig = scatter!([dpfe[1]], [mufe[1]], label = "CVM")
    return fig
end

function alocar(Sigma, mu, lista, ncarteiras = 10)
    wk = fef(Sigma, mu, ncarteiras)[1]
    ticklabel = "P" .* string.(collect(1:ncarteiras))
    ticklabel[1] = "CVM"
    ticklabel[ncarteiras] = "CRM"
    fig = groupedbar(wk', bar_position = :stack, bar_width=1.0, xlabel = "Carteiras", xticks=(1:ncarteiras, ticklabel), ylabel = "Pesos", label = permutedims(String.(lista)), legend = :outertopright, legendfontsize = 8, framestyle = :box)
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
