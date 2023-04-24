using PrecompileTools 

@setup_workload begin
    @compile_workload begin

      N = 8
      T = ComplexF32
      
      numCoils = 2
      I = T.(shepp_logan(N))
      I = circularShutterFreq!(I,1)
      
      coilsens = T.(birdcageSensitivity(N, numCoils, 1.5))
      
      # simulation parameters
      params = Dict{Symbol, Any}()
      params[:simulation] = "fast"
      params[:trajName] = "Spiral"
      params[:numProfiles] = 1
      params[:numSamplingPerProfile] = N*N
      params[:windings] = div(N,2)
      params[:AQ] = 3.0e-2
      params[:senseMaps] = coilsens
      
      # do simulation
      acqData = simulation(I, params)
      
      # reco parameters
      params = Dict{Symbol, Any}()
      params[:reco] = "multiCoil" #"standard"
      params[:reconSize] = (N,N)
      params[:regularization] = "L2"
      params[:iterations] = 2
      params[:solver] = "cgnr"
      params[:senseMaps] = coilsens
      

      reconstruction(acqData, params)
    end
end
