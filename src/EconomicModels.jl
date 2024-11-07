module EconomicModels

using Parameters, LinearSolve, SciMLOperators, Interpolations, ITensors, StaticArrays, Underscores
using EconomicSolvers
import Interpolations: coefficients
export BewleyAiyagariParameters, bewley_aiyagari_model

include("example_bewley_aiyagari.jl")

end
