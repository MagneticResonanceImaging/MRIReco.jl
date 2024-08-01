function MRIOperators.produ_inner!(K, C::matT, A::matT, shape, d::Vector{vecT}, s::vecT, sp, plan, idx, x_::vecT, p::Vector{arrT}) where {T, vecT <: AbstractGPUVector{T}, matT <: AbstractGPUMatrix{T}, arrT <: AbstractGPUArray{T}}
  for κ=1:K
    if !isempty(idx[κ])
      p[κ][:] .= C[κ,:] .* x_
      mul!(d[κ], plan[κ], p[κ])  
      # Assumption: l is unique per idx[κ]
      gpu_call(idx[κ], d[κ], view(A, :, κ), s) do ctx, indices, d_, A_, s_
        k = @linearidx(indices)
        l = indices[k]
        s_[l] += d_[k]*A_[l]
        return nothing
      end
    end
  end
  
  return
end


function MRIOperators.ctprodu_inner!(K, C::matT, A::matT, shape, d::Vector{vecT}, y::vecT, sp, plan, idx, x::vecT, p::Vector{arrT}) where {T, vecT <: AbstractGPUVector{T}, matT <: AbstractGPUMatrix{T}, arrT <: AbstractGPUArray{T}}

  for κ=1:K
    if !isempty(idx[κ])
      gpu_call(idx[κ], d[κ], view(A, :, κ), x) do ctx, indices, d_, A_, x_
        k = @linearidx(indices)
        l = indices[k]
        d_[k] = conj(A_[l]) * x_[l]
        return nothing
      end
      mul!(p[κ], adjoint(plan[κ]), d[κ])

      gpu_call(p[κ], view(C, κ, :), y) do ctx, p_, C_, y_
        k = @linearidx(p_)
        y_[k] += conj(C_[k]) * p_[k]
        return nothing
      end
    end
  end
    
  return
end