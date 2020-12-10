using PyPlot, MRIReco, FFTW, LinearAlgebra, BenchmarkTools

FFTW.set_num_threads(1);BLAS.set_num_threads(1)


threads = [1,2,4,8]
times = zeros(length(threads))

for t=1:length(times)
  @show threads[t]
  times[t] = parse(Float64,read(`julia-1.5 --threads=$(threads[t]) reconstruction.jl`, String))
end

plot(threads,times)
