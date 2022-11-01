module CompositeDyna

using CSV
using DelimitedFiles, DataFrames, LinearAlgebra
using Plots

#    include("CreateModel.jl")
include("FiniteElementModel.jl")
include("ForceStrainStiffnessMat.jl")
include("MassMatrix.jl")
include("ColumnHeight.jl")
include("Sky2Mat.jl")
include("BStiff.jl")

end # module CompositeDyna
