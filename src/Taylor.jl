struct Taylor{O,Nx,Ty}
    x::SVector{Nx,𝕣}
    y::Ty
end
Taylor{O}(x::SVector{Nx,𝕣},y::Ty) where{O,Nx,Ty<:AbstractArray} =       Taylor{O,Nx,Ty        }(x,y )
Taylor{O}(x::SVector{Nx,𝕣},y::Ty) where{O,Nx,Ty<:Tuple        } = Tuple(Taylor{O,Nx,typeof(yᵢ)}(x,yᵢ) for yᵢ∈y)
Taylor{O}(x::SVector{Nx,𝕣},y::Tuple{A    }) where{O,Nx,A    }   = tuple(Taylor{O}(x,y[1]))
Taylor{O}(x::SVector{Nx,𝕣},y::Tuple{A,B  }) where{O,Nx,A,B  }   = tuple(Taylor{O}(x,y[1]),Taylor{O}(x,y[2]))
Taylor{O}(x::SVector{Nx,𝕣},y::Tuple{A,B,C}) where{O,Nx,A,B,C}   = tuple(Taylor{O}(x,y[1]),Taylor{O}(x,y[2]),Taylor{O}(x,y[3]))
"""
    taylor = Taylor{O}(f,x₀)
    y      = taylor(x)

Compute the Taylor expansion of a function at `x₀` and evaluate it at `x`.

Input to `Taylor`:    
 - O ∈ {0,1,2}  is the order of the Taylor development.  It can be omitted

        taylor = Taylor(f,x₀)

    to compute a Taylor development of order equal to `min(2,precedence(x₀))`. 
 - `f` takes a single `SVector` as input.  Its output can
    be
    - a `Real`
    - a `SArray`
    - A `Tuple` which components are any combination of `Real`s and/or `SArray`s.  
      If for example `y,z = f(x)` then `taylor_y,taylor_z = Taylor{O}(f,x₀)`.  
 - `x₀` is the `SVector` at which the development is done.

`Taylor` objects can be called with a single argument: a `SVector` `x` of same length as `x₀`
at which to evaluate the Taylor expansion.

    y      = taylor(x)
    
"""
Taylor(f::Function,X::SVector{Nx,R}) where{Nx,R<:Real} = Taylor{min(precedence(R),2)}(f,X) 
function Taylor{0}(f::Function,X::SVector{Nx,R}) where{Nx,R<:Real}
    x  = VALUE(X)
    y  = f(x)
    return Taylor{0}(x,y)
end    
function Taylor{1}(f::Function,X::SVector{Nx,R}) where{Nx,R<:Real}
    x  = VALUE(X)
    y  = f(variate{1,Nx}(x))
    return Taylor{1}(x,y)
end    
function Taylor{2}(f::Function,X::SVector{Nx,R}) where{Nx,R<:Real}
    x  = VALUE(X)
    y  = f(variate{2,Nx}(variate{1,Nx}(x)))
    return Taylor{2}(x,y)
end    

(te::Taylor{O,Nx,Ty})(X::SVector{Nx,R}) where{O ,R<:Real,Nx,S,Ty<:SArray{S,𝕣}} = te.y  # y[] was untouched by x or to 0th order 
(te::Taylor{O,Nx,Ty})(X::SVector{Nx,R}) where{O ,R<:Real,Nx,  Ty<:         𝕣 } = te.y  # y   was untouched by x or to 0th order 
(te::Taylor{0,Nx,Ty})(X::SVector{Nx,R}) where{Nx,R<:Real     ,Ty             } = te.y  

function (te::Taylor{1,Nx,Ty})(X::SVector{Nx,R}) where{R<:Real,dR<:∂ℝ,S,Nx,Ty<:SArray{S,dR}} # y[] to 1st order
    ΔX = X-te.x
    SArray{S}(yᵢ.x + yᵢ.dx∘₁ΔX  for yᵢ∈ te.y)
end
function (te::Taylor{1,Nx,Ty})(X::SVector{Nx,R}) where{R<:Real,Nx,Ty<:∂ℝ}                  # y  to 1st order
    ΔX = X-te.x
    te.y.x + te.y.dx∘₁ΔX
end

function (te::Taylor{2,Nx,Ty})(X::SVector{Nx,R}) where{R<:Real,dR<:∂ℝ,S,Nx,Ty<:SArray{S,dR}}          # y[] to 2nd order
    ΔX = X-te.x   
    SArray{S}(yᵢ.x.x + sum(    (yᵢ.dx[j].x + .5*(yᵢ.dx[j].dx ∘₁ ΔX)) *ΔX[j] for j∈eachindex(ΔX))  for yᵢ∈ te.y)  
end
function (te::Taylor{2,Nx,Ty})(X::SVector{Nx,R}) where{R<:Real,Nx,Ty<:∂ℝ}                 # y  to 2nd order 
    ΔX = X-te.x   
    te.y.x.x + sum(    (te.y.dx[j].x + .5*(te.y.dx[j].dx ∘₁ ΔX)) *ΔX[j] for j∈eachindex(ΔX))   
end

∂(t::Taylor{O,Nx}) where{O,Nx} = Taylor{O-1}(t.x,∂{O,Nx}(t.y))
∂(t::Taylor{0   })             = muscadeerror("Tried to differentiate a 0th order Taylor expansion")

# TODO: return a tuple of expansions, not an expansion of a tuple  DONE
# TODO: Taylor stores Adiff objects. 
# TODO Dots: no need to operate on Tuples DONE

