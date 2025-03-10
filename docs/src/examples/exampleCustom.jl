using MRIReco, MRIFiles, MRIOperators, MRISampling
using MRIReco.MRIOperators.Wavelets, MRIReco.RegularizedLeastSquares
using DelimitedFiles, LinearAlgebra, CairoMakie
include(joinpath(@__DIR__,"exampleUtils.jl"))

function analyzeImage(res, x::Vector{T},D::Matrix{T},xsize::NTuple{2,Int64},psize::NTuple{2,Int64};t0::Int64=size(D,2),tol=1e-3) where T
  nx,ny = xsize
  px,py = psize
  x = reshape(x,nx,ny)
  x_pad = repeat(x,2,2)[1:nx+px-1,1:ny+py-1] # pad image using periodic boundary conditions
  α = zeros(T,size(D,2),nx,ny)
  patch = zeros(T,px*py)
  for j=1:ny
    for i=1:nx
      patch[:] .= vec(x_pad[i:i+px-1,j:j+py-1])
      norm(patch)==0 && continue
      # matchingpursuit is contained in Wavelets.jl
      α[:,i,j] .= matchingpursuit(patch, x->D*x, x->transpose(D)*x, tol)
    end
  end
  return res .= vec(α)
end

function synthesizeImage(res, α::Vector{T},D::Matrix{T},xsize::NTuple{2,Int64},psize::NTuple{2,Int64}) where T
  nx,ny = xsize
  px,py = psize
  x = zeros(T,nx+px-1,ny+py-1)
  α = reshape(α,:,nx,ny)
  for j=1:ny
    for i=1:nx
      x[i:i+px-1,j:j+py-1] .+= reshape(D*α[:,i,j],px,py)
    end
  end
  return res .= vec(x[1:nx,1:ny])/(px*py)
end

function dictOp(D::Matrix{T},xsize::NTuple{2,Int64},psize::NTuple{2,Int64},tol::Float64=1.e-3) where T
        produ = (res,x)->analyzeImage(x,D,xsize,psize,tol=tol)
        ctprodu = (res,x)->synthesizeImage(x,D,xsize,psize)
  return LinearOperator(T,prod(xsize)*size(D,2),prod(xsize),false,false
          , produ
          , nothing
          , ctprodu )
end
## To test our method, let us load some simulated data and subsample it

# phantom
img = readdlm(joinpath(@__DIR__,"data/mribrain100.tsv"))

acqData = AcquisitionData(ISMRMRDFile(joinpath(@__DIR__,"data/acqDataBrainSim100.h5")))
nx,ny = acqData.encodingSize

# undersample kspace data
acqData = sample_kspace(acqData, 2.0, "poisson", calsize=25,profiles=false);


# CS reconstruction using Wavelets
params = Dict{Symbol,Any}()
params[:reco] = "direct"
params[:reconSize] = (nx,ny)

@time img_u = reconstruction(acqData,params)
@info "relative error: $(norm(img-img_u)/norm(img))"

## Load the Dictionary, build the sparsifying transform and perform the reconstruction

# load the dictionary
D = ComplexF32.(readdlm(joinpath(@__DIR__,"data/brainDict98.tsv")))

begin
D2 = reshape(D,6,6,:)# check patch
f=Figure()
for i in 1:6, j in 1:6
  ax=Axis(f[i,j],aspect=1)
  heatmap!(ax,abs.(D2[:,:,(i-1)*6+j]))
  hidedecorations!(ax)
end
f
end

# some parameters
px, py = (6,6)  # patch size
K = px*py       # number of atoms

# CS reconstruction using Wavelets
params = Dict{Symbol,Any}()
params[:reco] = "standard"
params[:reconSize] = (nx,ny)
params[:iterations] = 50
params[:reg] = L1Regularization(0.0)
params[:sparseTrafo] = dictOp(D,(nx,ny),(px,py),2e-2)
params[:rho] = 0.1
params[:solver] = ADMM
params[:absTol] = 1.e-4
params[:relTol] = 1.e-2
params[:normalizeReg] = MeasurementBasedNormalization()
params[:kwargWarning] = true

@time img_d = reconstruction(acqData,params)
@info "relative error: $(norm(img-img_d)/norm(img))"


# use CairoMakie for interactive display
begin
  colormap=:grays
  f = Figure(size=(1000,300))
  ax = Axis(f[1,1],title="Original")
  heatmap!(ax,rotr90(abs.(img[:,:,1,1,1]));colormap)
  ax = Axis(f[1,2],title="Undersampled")
  heatmap!(ax,rotr90(abs.(img_u[:,:,1,1,1]));colormap)
  ax = Axis(f[1,3],title="Wavelet")
  heatmap!(ax,rotr90(abs.(img_w[:,:,1,1,1]));colormap)
  ax = Axis(f[1,4],title="Dictionary")
  heatmap!(ax,rotr90(abs.(img_d[:,:,1,1,1]));colormap)
  [hidedecorations!(f.content[ax]) for ax in eachindex(f.content)]
  f
end
# For comparison, let us perform the same reconstruction as above but with a Wavelet transform

params = Dict{Symbol,Any}()
params[:reco] = "standard"
params[:reconSize] = (nx,ny)
params[:iterations] = 50
params[:reg] = L1Regularization(2.e-2)
params[:sparseTrafo] = "Wavelet"

img_w = reconstruction(acqData,params)
@info "relative error: $(norm(img-img_w)/norm(img))"


# use CairoMakie for interactive display
begin
  colormap=:grays
  f = Figure(size=(1000,300))
  ax = Axis(f[1,1],title="Original")
  heatmap!(ax,rotr90(abs.(img[:,:,1,1,1]));colormap)
  ax = Axis(f[1,2],title="Undersampled")
  heatmap!(ax,rotr90(abs.(img_u[:,:,1,1,1]));colormap)
  ax = Axis(f[1,3],title="Wavelet")
  heatmap!(ax,rotr90(abs.(img_w[:,:,1,1,1]));colormap)
  ax = Axis(f[1,4],title="Dictionary")
  heatmap!(ax,rotr90(abs.(img_d[:,:,1,1,1]));colormap)
  [hidedecorations!(f.content[ax]) for ax in eachindex(f.content)]
  f
end

# export images
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/brainOrig.png")
exportImage(filename, abs.(img[:,:,1,1,1]) )
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/brainWavelet.png")
exportImage(filename, abs.(img_w[:,:,1,1,1]) )
filename = joinpath(dirname(pathof(MRIReco)),"../docs/src/assets/brainDict.png")
exportImage(filename, abs.(img_d[:,:,1,1,1]) )
