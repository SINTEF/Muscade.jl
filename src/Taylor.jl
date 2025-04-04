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



##############


struct motion{P}          end 
struct motion_{P,Q}       end 
struct motion⁻¹{P,ND,OD}  end 
"""
    P  = constants(X,U,A,t)
    X_ = motion{P}(X)

Transform a `NTuple` of `SVector`s, for example the vector `X` provided as an input to
`residual` or `Lagrangian` into a `SVector` of `∂ℝ`.  This can be used by an element to 
compute time derivatives, for example Euler, Coriolis and centrifugal accelerations, 
or strain rates.

Some principles of safe automatic differentiation must be adhered to:
- the function that uses `motion` must also 'unpack' : no variable that is touched by 
  the output of `motion` must be returned by the function without having been unpacked
  by `motion⁻¹`.
- The precendence `P` must be calculated using `constants` with all variables that are input to 
  the function and may be differentiated.
- If other levels of automatic differentiation are introduced within the function, unpack in reverse
  order of packing.    

See [`motion⁻¹`](@ref)
"""
motion{ P    }(a::NTuple{ND,SV{N,R}}) where{ND,P,N,R        } = SV{N}(motion_{P,P+ND-1}(ntuple(j->a[j][i],ND)) for i=1:N)
motion_{P,Q  }(a::NTuple{D,      R }) where{D ,P  ,R<:Real,Q} = ∂ℝ{Q,1}(motion_{P,Q-1}(a),SV(motion_{P,Q-1}(a[2:D]))) 
motion_{P,P  }(a::NTuple{D,      R }) where{D ,P  ,R<:Real  } = a[1]

"""
    P  = constants(X,U,A,t)
    ND = length(X)
    X_  = motion{P,ND}(X)
    Y_  = f(Y_)    
    Y₀ = motion⁻¹{P,ND,0}(Y_)
    Y₁ = motion⁻¹{P,ND,1}(Y_)
    Y₂ = motion⁻¹{P,ND,2}(Y_)
    Y  = motion⁻¹{P,ND  }(Y) 

Extract the value and time derivatives from a variable that is a function of the output of `motion`.
In the above `Y` is a tuple of same length as `X`.  One can use `∂0`,`∂1` and `∂2` to unpack `Y`.    

See [`motion`](@ref)
"""
motion⁻¹{P,1,0}(a::ℝ) where{P} =                             a
motion⁻¹{P,2,0}(a::ℝ) where{P} =              value{P+2-1  }(a)
motion⁻¹{P,3,0}(a::ℝ) where{P} = value{P+3-2}(value{P+3-1  }(a))
# velocities
motion⁻¹{P,1,1}(a::ℝ) where{P} = 0. 
motion⁻¹{P,2,1}(a::ℝ) where{P} =              ∂{    P+2-1,1}(a)[1]  # [1]: only partial is wrt time
motion⁻¹{P,3,1}(a::ℝ) where{P} = value{P+3-2}(∂{    P+3-1,1}(a)[1])
# accelerations
motion⁻¹{P,1,2}(a::ℝ) where{P} = 0. 
motion⁻¹{P,2,2}(a::ℝ) where{P} = 0.
motion⁻¹{P,3,2}(a::ℝ) where{P} = ∂{   P+3-2,1}(∂{   P+3-1,1}(a)[1])[1]

motion⁻¹{P,ND,OD}(a::AbstractArray) where{P,ND,OD} = motion⁻¹{P,ND,OD}.(a)
#motion⁻¹{P,ND   }(a               ) where{P,ND   } = ntuple(ID->motion⁻¹{P,ND,ID-1}(a) ,ND)
motion⁻¹{P,1    }(a               ) where{P   } = (motion⁻¹{P,1,0}(a),)
motion⁻¹{P,2    }(a               ) where{P   } = (motion⁻¹{P,2,0}(a),motion⁻¹{P,2,1}(a))
motion⁻¹{P,3    }(a               ) where{P   } = (motion⁻¹{P,3,0}(a),motion⁻¹{P,3,1}(a),motion⁻¹{P,3,2}(a))



