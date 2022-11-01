using CSV
using DataFrames
#using Plots
#using BendingStiff
using DelimitedFiles
#using LinearAlgebra
#using CompositeDyna

data = readdlm("./test/data.in", header=false)
df = DataFrame(data, :auto)

iatype::Int64 = df[1, "x1"]
nlayers::Int64 = df[2, "x1"]
E11::Float64 = df[2, "x2"]
E22::Float64 = df[2, "x3"]
PR12::Float64 = df[2, "x4"]
G12::Float64 = df[2, "x5"]
G13::Float64 = df[2, "x6"]
G23::Float64 = df[2, "x7"]

ori = Vector{Float64}(undef, nlayers)
zpos = Vector{Float64}(undef, nlayers+1)
for i=1:nlayers
    ori[i] = df[2+i, "x1"]
    zpos[i] = df[2+nlayers+i, "x1"]
end
zpos[nlayers+1] = df[2+2*nlayers+1, "x1"]

th::Float64=zpos[1]-zpos[nlayers+1]

# Length - ac, Width - bc, No. of elements along length - nx
# No. of elements along width - ny
ac::Float64=df[2+2*nlayers+2, "x1"]
bc::Float64=df[2+2*nlayers+2, "x2"]
nx::Int64=df[2+2*nlayers+2, "x3"]
ny::Int64=df[2+2*nlayers+2, "x4"]

# Read actual boundary conditions on left, right, top and bottom
kl::Int64=df[2+2*nlayers+3, "x1"]
kr::Int64=df[2+2*nlayers+3, "x2"]
kt::Int64=df[2+2*nlayers+3, "x3"]
kb::Int64=df[2+2*nlayers+3, "x4"]

# Read plane stress boundary conditions on left, right, top and bottom
klp::Int64=df[2+2*nlayers+4, "x1"]
krp::Int64=df[2+2*nlayers+4, "x2"]
ktp::Int64=df[2+2*nlayers+4, "x3"]
kbp::Int64=df[2+2*nlayers+4, "x4"]