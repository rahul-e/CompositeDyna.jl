using CSV
using DataFrames
#using Plots
#using BendingStiff
using Libdl
using DelimitedFiles
using LinearAlgebra
using CompositeDyna

# Load the Fortran library
const libsubspace= dlopen("./include/libsubspace.so")

libfunc = dlsym(libsubspace, :subspace_)


# Define the ccall to interface with the Fortran subroutine
function subspace(
    gk::Vector{Float64}, gm::Vector{Float64}, nds::Vector{Int64},
    r::Array{Float64,2}, eigv::Vector{Float64}, tt::Vector{Float64}, 
    w::Vector{Float64}, vec::Array{Float64,2}, d::Vector{Float64}, 
    rtolv::Vector{Float64}, nn::Int64, nnm::Int64, 
    nwk::Int64, nwm::Int64, nroot::Int64, 
    rtol::Float64, nc::Int64, nnc::Float64,
    nitem::Int64, ifss::Int64, ifpr::Int64, 
    nrukku::Int64, ddd::Array{Float64,2})
    ccall(libfunc,
          Cvoid,                 # Return type
          (Ref{Float64}, Ref{Float64}, Ref{Int64}, 
           Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Float64}, Ref{Float64},
           Ref{Float64}, Ref{Int64}, Ref{Int64},
           Ref{Int64}, Ref{Int64}, Ref{Int64},
           Ref{Float64}, Ref{Int64}, Ref{Float64},
           Ref{Int64}, Ref{Int64}, Ref{Int64},
           Ref{Int64}, Ref{Float64}),
           gk, gm, nds, 
           r, eigv, tt,
           w, vec, d,
           rtolv, nn, nnm,
           nwk, nwm, nroot,
           rtol, nc, nnc, 
           nitem, ifss, ifpr,
           nrukku, ddd)  # Arguments passed to the subroutine
    return
end


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

ori = fill(0.0, nlayers)
zpos = fill(0.0, nlayers+1)
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

println("Total number of elements ",nx*ny)

CC = CompositeDyna.composit(nlayers, ori, zpos, E11, E22, PR12, G12, G13, G23)

gk, gm, nbig, nsky, cht, nds = CompositeDyna.geo(CC, ac, bc, th, nx, ny, kl, kr, kt, kb, klp, krp, ktp, kbp)

println("NDS[end], CHT[end] = ", nds[end]," ", cht[end])

GK = CompositeDyna.Sky2Mat(gk, nbig, nsky, cht)
GM = CompositeDyna.Sky2Mat(gm, nbig, nsky, cht)

println("NBIG ", nbig)
# Specific Weight
spwt::Float64 = 2.8E-5
# Density
rho::Float64 = spwt/9810

if (iatype==1) 
    freqn=bc*bc*sqrt(rho/(E11*th^3/(12*(1-PR12^2))))
elseif (iatype==2)
    freqn=bc*bc*sqrt(rho*th/(E22*th^2))
end

open("GK_mat.txt", "w") do file
    for i in 1:nbig
         for j in 1:nbig
         print(file, GK[i,j], " ")
         end
         println(file, " ")
    end
end

# Parameters controlling the Eigen solver 
nrukku = 0
nn = nbig
nnm = nbig + 1
nwk = nsky
nwm = nwk
nroot = 40
rtol = 0.00001
nc = 48
nnc = nc*(nc+1)/2
nitem = 16
ifss = 1
ifpr = 1
r = fill(0.0, (nn, nc))
eigv = fill(1.0, nc)
tt = fill(0.0, nn)
w = fill(0.0, nn)
vec = fill(0.0, (nc, nc))
d = fill(0.0, nc)
rtolv = fill(0.0, nc)
ddd = fill(0.0, (nc, 6000))

#subspace(nrukku, nnm, nwm, nroot, rtol, nc, nnc, nitem, ifss, ifpr, nbig, nsky)
#subspace(gk, gm, nds, r, eigv, tt, w, vec, d, rtolv, nn, nnm, nwk, nwm, nroot, rtol, nc, 
#nnc, nitem, ifss, ifpr, nrukku, ddd)
values = eigvals(GK, GM)
for i=1:2
    println(values[i])
end

