function MRIOperators.produ_inner!(K, C::matT, A::matT, shape, d::Vector{vecT}, s::vecT, sp, plan, idx, x_::vecT, p::Vector{arrT}) where {T, vecT <: AbstractGPUVector{T}, matT <: AbstractGPUMatrix{T}, arrT <: AbstractGPUArray{T}}
  
  @kernel inbounds = true cpu = false function produ_inner_kernel!(indices, d, A, s)
    k = @index(Global, Linear)
    # Assumption: l is unique per idx[κ]
    l = indices[k]
    s[l] += d[k]*A[l]
  end
  kernel! = produ_inner_kernel!(get_backend(C))

  for κ=1:K
    if !isempty(idx[κ])
      p[κ][:] .= C[κ,:] .* x_
      mul!(d[κ], plan[κ], p[κ])
      kernel!(idx[κ], d[κ], view(A, :, κ), s; ndrange = length(idx[κ]))
    end
  end
  
  return
end


function MRIOperators.ctprodu_inner!(K, C::matT, A::matT, shape, d::Vector{vecT}, y::vecT, sp, plan, idx, x::vecT, p::Vector{arrT}) where {T, vecT <: AbstractGPUVector{T}, matT <: AbstractGPUMatrix{T}, arrT <: AbstractGPUArray{T}}

  
  @kernel inbounds = true cpu = false function ctprodu_inner_1!(indices, d, A, x)
    k = @index(Global, Linear)
    l = indices[k]
    d[k] = conj(A[l]) * x[l]
  end
  kernel1! = ctprodu_inner_1!(get_backend(C))

  @kernel inbounds = true cpu = false function ctprodu_inner_2!(p, C, y)
    k = @index(Global, Linear)
    y[k] += conj(C[k]) * p[k]
  end
  kernel2! = ctprodu_inner_2!(get_backend(C))


  for κ=1:K
    if !isempty(idx[κ])
      kernel1!(idx[κ], d[κ], view(A, :, κ), x; ndrange = length(idx[κ]))
      mul!(p[κ], adjoint(plan[κ]), d[κ])
      kernel2!(p[κ], view(C, κ, :), y; ndrange = length(p[κ]))
    end
  end
    
  return
end