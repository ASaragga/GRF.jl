module GRF

using ARCHModels                   
using CSV
using Dates             
using DataFrames
using Distributions 
using GLM  
using HypothesisTests   
using HTTP
using Ipopt             
using JSON
using JuMP 
using LaTeXStrings                       
using LinearAlgebra  
using LsqFit  
using Random           
using Statistics        
using StatsPlots

export yahoo, retornos, Dates
export SP500n

include("modulo00.jl")
include("modulo02.jl")
include("modulo03.jl")

# Notebooks interativos 
# binder: https://mybinder.org 
# https://mybinder.org/v2/gh/ASaragga/GRF.jl/HEAD
#           
# nbviewer: https://nbviewer.jupyter.org
# https://github.com/ASaragga/GRF.jl/blob/master/src/Modulo1.1-1.7.ipynb


end # module
