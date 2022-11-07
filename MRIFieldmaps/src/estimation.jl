"""
ml_cost(x::Matrix{T},y::Matrix{Complex{T}},z::Matrix{Complex{T}}, β)
Calculates the regularized ML estimator cost between an estimated phase map and multi-echo scan data
# Arguments
* `x::Matrix{T}` - fieldmap estimate (in radians)
* `y::Matrix{Complex{T}}` - Complex first-echo image data
* `z::Matrix{Complex{T}}` - Complex second-echo image data
* `m::Matrix{T}` - normalized weighting data (m ∈ [0,1]:= abs.(conj.(y).*z)./maximum(abs.(conj.(y).*z))), precomputed for speed
* `β` - Regularization parameter controlling roughness penalty
"""
function ml_cost(x::Matrix{T},y::Matrix{Complex{T}},z::Matrix{Complex{T}},m::Matrix{T}, β) where T

    Ψ = (sum(m.* (1 .- cos.(angle.(z) .- angle.(y) .- x)),dims=[1,2]) + β * R(x))[1]

end

"""
R(x::Matrix{T})
Regularization term which penalizes roughness in a first-order finite-difference sense
# Arguments
* `x::Matrix{T}` - fieldmap estimate (in radians)
"""
function R(x::Matrix{T}) where T

    return T.(0.5) * sum(abs2,diff(x,dims=1),dims=[1,2]) + T.(0.5) * sum(abs2,diff(x,dims=2),dims=[1,2])

end

"""
pcg_ml_est_fieldmap(y::AbstractMatrix{Complex{T}},z::AbstractMatrix{Complex{T}},β) 
Estimates the fieldmap (low-level) using the method presented in https://doi.org/10.1109/tmi.2008.923956
# Required Arguments
* `y::AbstractMatrix{Complex{T}}` - Complex first-echo image data
* `z::AbstractMatrix{Complex{T}}` - Complex second-echo image data
# Optional Arguments
* `β` - Regularization parameter controlling roughness penalty
* `reltol` - early stopping criteria (exit if subsequent cost function change < reltol)
"""
function pcg_ml_est_fieldmap(y::AbstractMatrix{Complex{T}},z::AbstractMatrix{Complex{T}},β = 1e-3, reltol = 5e-3; quiet = true) where T

    # precomputation of static cost function inputs
    d = 1
    κ = conj.(y) .* z
    x = angle.(κ)
    m = abs.(κ) ./ maximum(abs.(κ))

    # ideally this is SNR dependent trust, not signal dependent trust, but signal is a reasonable proxy for SNR
    trust_step = 1.0 ./ (m .+ 4*β)

    c = 0
    Δ = 1
    itcount = 0

    # gradient descent routine
    while Δ > reltol

        gs = gradient(Flux.params(x)) do
            c = ml_cost(x,y,z,m,β) # TODO: need to interpolate the y z and β to have better performance 
        end

        x̄ = gs[x]
        x .-= (trust_step) .* x̄ 

        Δ = abs.(ml_cost(x,y,z,m,β) - c)/c

        itcount +=1

    end

    if !quiet
        print("Required $itcount iterations to converge below tolerance $reltol ")
    end

    return x

end

"""
estimate_cmap(imData,slices, TE1,TE2,β,isrotated)
Processes 5D volume data (high-level) as output from MRIReco.reconstruction to estimate fieldmaps using the method presented by Funai and Fessler
# Required Arguments
* `imData` - 5-D array with complex image data -> [X,Y,Z,ECHO,REPETITION]
* `slices` - vector of slices to process (must be within range of 3rd dimension of imData)
* `TE1` - Echo time 1 [ms]
* `TE2` - Echo time 2 [ms]
# Optional Arguments
* `isrotated` - Boolean controlling whether to rotate the B0 maps to match the images or not (legacy feature)
# Keyword Arguments 
* `β` - Regularization parameter controlling roughness penalty, appropriate range: β ∈ [1e-5,1.0]
* `reltol` - early stopping criteria (exit if ((ϕₙ₊₁ - ϕₙ)/ϕₙ) < reltol)
"""
function estimate_cmap(imData,slices, TE1,TE2,isrotated; β = 1e-3, reltol = 0.001)

    cmap = Complex.(zeros(size(imData)[1:3]))

    @sync for slice in slices
        Threads.@spawn begin
            cmap[:,:,slice] = pcg_ml_est_fieldmap(imData[:,:,slice,1,1],imData[:,:,slice,2,1],β,reltol) ./ ((TE2 - TE1)/1000)
            #println("for slice $slice")
        end
    end

    cmap = mapslices(isrotated ? x->x : x-> rotl90(x),cmap,dims=(1,2)) 

    return real.(cmap) # can return complex but not supported yet

end