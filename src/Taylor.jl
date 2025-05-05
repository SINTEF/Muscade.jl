struct motion{P}          end 
struct motion_{P,Q}       end 
struct motion⁻¹{P,ND,OD}  end 
"""
    P  = constant(X)
    X_ = motion{P}(X)

Transform a `NTuple` of `SVector`s, for example the vector `X` provided as an input to
`residual` or `Lagrangian` into a `SVector` of `∂ℝ`.  This can be used by an element to 
compute time derivatives, for example Euler, Coriolis and centrifugal accelerations, 
or strain rates.

Some principles of safe automatic differentiation must be adhered to:
- the function that uses `motion` must also 'unpack' : no variable that is touched by 
  the output of `motion` must be returned by the function without having been unpacked
  by `motion⁻¹`. Touched variables can for example be marked with an underscore
- The precendence `P` must be calculated using `constants` with all variables that are input to 
  the function and may be differentiated.
- If other levels of automatic differentiation are introduced within the function, unpack in reverse
  order of packing.    

See [`motion⁻¹`](@ref)
"""
motion{ P    }(a::NTuple{ND,SV{N,R}}) where{ND,P,N,R        } = SV{N}(motion_{P-1,P+ND-2}(ntuple(j->a[j][i],ND)) for i=1:N)
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
In the above `Y` is a tuple of length `ND`.  One can use `∂0`,`∂1` and `∂2` to unpack `Y`.    

See also [`motion`](@ref)
"""
motion⁻¹{P,1,0}(a::ℝ) where{P} =                             a
motion⁻¹{P,2,0}(a::ℝ) where{P} =            value{P+1  }(a)
motion⁻¹{P,3,0}(a::ℝ) where{P} = value{P+1}(value{P+2 }(a))
# velocities
motion⁻¹{P,1,1}(a::ℝ) where{P} = 0. 
motion⁻¹{P,2,1}(a::ℝ) where{P} =            ∂{    P+1  ,1}(a)[1]  # [1]: only partial is wrt time
motion⁻¹{P,3,1}(a::ℝ) where{P} = value{P+1}(∂{    P+2,1}(a)[1])
# accelerations
motion⁻¹{P,1,2}(a::ℝ) where{P} = 0. 
motion⁻¹{P,2,2}(a::ℝ) where{P} = 0.
motion⁻¹{P,3,2}(a::ℝ) where{P} = ∂{   P+1,1}(∂{   P+2,1}(a)[1])[1]

motion⁻¹{P,ND,OD}(a::AbstractArray) where{P,ND,OD} = motion⁻¹{P,ND,OD}.(a)
#motion⁻¹{P,ND   }(a               ) where{P,ND   } = ntuple(ID->motion⁻¹{P,ND,ID-1}(a) ,ND)
motion⁻¹{P,1    }(a               ) where{P   } = (motion⁻¹{P,1,0}(a),)
motion⁻¹{P,2    }(a               ) where{P   } = (motion⁻¹{P,2,0}(a),motion⁻¹{P,2,1}(a))
motion⁻¹{P,3    }(a               ) where{P   } = (motion⁻¹{P,3,0}(a),motion⁻¹{P,3,1}(a),motion⁻¹{P,3,2}(a))


#############

struct revariate{ O,N}   end
struct revariate_{P,N}   end
revariate(a) = revariate{0}(a)
""" 
    TX = revariate{O}(X)

The vector `X` of `Real`s (possibly: `∂ℝ`s) is stripped of its partials, an revariated to the
order `precedence(X)+O`.

    TX = revariate(X)

revariates to the order `precedence(X)`.  

`revariate`, in conjunction with `compose` can be used to improve performance when the length of 
`X` is smaller than the length of its partials.

Be extremely careful never to mix any variable that is a function of `X` with any other variables
containing  `∂ℝ`s but not produced by the same `revariate`.

See also: [`compose`](@ref)
"""
function revariate{O}(a::SV{N,R}) where{O,N,R} 
    P  = precedence(R)+O
    va = VALUE(a)
    P==0 ? va : SV{N}(∂ℝ{P,N  }(revariate_{P-1,N}(va[i],i),i) for i=1:N)
end
revariate_{P,N}(a,i) where{P,N} = ∂ℝ{P,N  }(revariate_{P-1,N}(a,i),i)
revariate_{0,N}(a,i) where{  N} =                             a

"""
    compose(Ty,Δx)

`Ty::∂ℝ` has partials to arbitrary order with respect to a variable `x`. These
partials define a Taylor expansion, which `compose` evaluates at an increment `Δx` 
of the value `x` at which `Ty` was computed.

`compose` also handles nested structures of `Tuple`s and `SVector`s of `∂ℝ`, applying the
composition to each element.

`compose` can be used to evaluate a Taylor expansion, but it can also be used to compose
automatic differentiation.  For example 
    Tx = revariate(x)
    Ty = f(Tx)
    y  = compose(Ty,x-VALUE(x))    
is faster than
    y  = f(x)
if the length of `x` is smaller than the length of its partials.

See also: [`revariate`](@ref), [`fast`](@ref)    
"""
compose(y::Tuple,Δx)                          = tuple(compose(first(y),Δx),compose(Base.tail(y),Δx)...) 
compose(y::Tuple{},Δx)                        = tuple() 
compose(y::SArray{S},Δx) where{S}             = SArray{S}(compose(yᵢ,Δx) for yᵢ∈y) 
compose(y::∂ℝ,Δx)                             = compose(y.x,Δx) + compose_right(y,Δx)
compose(y::𝕣 ,Δx)                             =         y
compose_right(y::∂ℝ{P},Δx::SVector{N}) where{P,N} = sum(compose_right(y.dx[i],Δx)*Δx[i] for i∈1:N)*(1/P)
compose_right(y::𝕣    ,Δx            )            =               y

"""
    y,... = fast(f,x)

In the context of forward automatic differentiation using `∂ℝ`, accelerate the evaluation of
`y,...= f(x)` if the length of `x` is smaller than the length of its partials.

Be extremely careful with closures, making sure that `f` does not capture variables of type `∂ℝ`.

Wrapper function of [`revariate`](@ref) and [`compose`](@ref)      
"""
fast(f,x) = compose(f(revariate(x)),x-VALUE(x))    

"""
    composewithJacobian{P,ND,NDOF}

Works, but still work to do on the syntactic sugar.    
"""
struct composewithJacobian{P,ND,NDOF} end
function composewithJacobian{P,ND,NDOF}(Ty,X_) where{P,ND,NDOF}
    X₀         = motion⁻¹{P-1,ND,0}(X_)
    y          = motion⁻¹{P-1,ND  }(compose(value{P}( Ty  ),X_))
    y∂X₀       =                  compose(∂{P,NDOF}(Ty  ),X₀ )
    return y,y∂X₀
end
