
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
    y      = taylor(x₁)

or    

    y      = Taylor{O}(f,x₀)(x₁)      

 - O ∈ {0,1,2}  is the order of the Taylor development
 - `f` must be a `SVector`-valued function of a  `SVector`
 - `x₀` is the `SVector` at which the development is done.
 
    y      = Taylor(f,x₀)(x₁) 
    
(without specifying the order) computes a Taylor development of order equal to `precedence(x₀)`    
 
"""
# TODO: return a tuple of expansions, not an expansion of a tuple  DONE
# TODO: Taylor stores Adiff objects. 
# TODO Dots: no need to operate on Tuples DONE

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

(te::Taylor{0,Nx,Ty})(X::SVector{Nx,R}) where{Nx,R<:Real,Ty} = te.y 
function (te::Taylor{1,Nx,Ty})(X::SVector{Nx,R}) where{R<:Real,S,Nx,Ty<:SArray{S}}
    ΔX = X-te.x
    SArray{S}(yᵢ.x + yᵢ.dx∘₁ΔX  for yᵢ∈ te.y)
end
function (te::Taylor{2,Nx,Ty})(X::SVector{Nx,R}) where{R<:Real,S,Nx,Ty<:SArray{S}}
    ΔX = X-te.x   
    SArray{S}(yᵢ.x.x + sum(    (yᵢ.dx[j].x + .5*(yᵢ.dx[j].dx ∘₁ ΔX)) *ΔX[j] for j∈eachindex(ΔX))  for yᵢ∈ te.y)  
end

∂(t::Taylor{O,Nx}) where{O,Nx} = Taylor{O-1}(t.x,∂{O,Nx}(t.y))
∂(t::Taylor{0   }) where{O   } = MuscadeError("Tried to differentiate a 0th order Taylor")

