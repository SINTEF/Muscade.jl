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
motion⁻¹{P,2,0}(a::ℝ) where{P} =          value{P   }(a)
motion⁻¹{P,3,0}(a::ℝ) where{P} = value{P}(value{P+1 }(a))
# velocities
motion⁻¹{P,1,1}(a::ℝ) where{P} = 0. 
motion⁻¹{P,2,1}(a::ℝ) where{P} =          ∂{    P  ,1}(a)[1]  # [1]: only partial is wrt time
motion⁻¹{P,3,1}(a::ℝ) where{P} = value{P}(∂{    P+1,1}(a)[1])
# accelerations
motion⁻¹{P,1,2}(a::ℝ) where{P} = 0. 
motion⁻¹{P,2,2}(a::ℝ) where{P} = 0.
motion⁻¹{P,3,2}(a::ℝ) where{P} = ∂{   P,1}(∂{   P+1,1}(a)[1])[1]
motion⁻¹{P,ND,OD}(a::SArray{S,R}) where{S,P,ND,OD,R<:ℝ}   = SArray{S}(motion⁻¹{P,ND,OD}(aᵢ) for aᵢ∈a)

motion⁻¹{P,1    }(a::Union{ℝ,SArray}) where{P   } = (motion⁻¹{P,1,0}(a),)
motion⁻¹{P,2    }(a::Union{ℝ,SArray}) where{P   } = (motion⁻¹{P,2,0}(a),motion⁻¹{P,2,1}(a))
motion⁻¹{P,3    }(a::Union{ℝ,SArray}) where{P   } = (motion⁻¹{P,3,0}(a),motion⁻¹{P,3,1}(a),motion⁻¹{P,3,2}(a))
motion⁻¹{P,ND   }(a::Union{Tuple,NamedTuple}) where{P,ND} = map(motion⁻¹{P,ND},a)
motion⁻¹{P,ND   }(a...)               where{P,ND} = motion⁻¹{P,ND}(a)

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
    McLaurin(Ty,x)

`Ty::∂ℝ` has partials to arbitrary order with respect to a variable `x`. These
partials define a McLaurin expansion, which `McLaurin` evaluates at value `x`, 
as if `Ty` had been computed at 0.

`McLaurin` handles nested structures of `Tuple`s and `SVector`s of `∂ℝ`, applying the
expansion to each element.

`McLaurin` is a utility function behind [`compose`](@ref) and [`Taylor`](@ref)

See also: [`compose`](@ref), [`Taylor`](@ref), [`revariate`](@ref), [`fast`](@ref)    
"""
McLaurin(y::Tuple,Δx)                          = tuple(McLaurin(first(y),Δx),McLaurin(Base.tail(y),Δx)...) 
McLaurin( ::Tuple{},Δx)                        = tuple() 
McLaurin(y::SArray{S},Δx) where{S}             = SArray{S}(McLaurin(yᵢ,Δx) for yᵢ∈y) 
McLaurin(y::∂ℝ,Δx)                             = McLaurin(y.x,Δx) + McLaurin_right(y,Δx)
McLaurin(y::𝕣 ,Δx)                             =          y
McLaurin_right(y::∂ℝ{P},Δx::SVector{N}) where{P,N} = sum(McLaurin_right(y.dx[i],Δx)*Δx[i] for i∈1:N)*(1/P)
McLaurin_right(y::𝕣    ,Δx            )        =          y

"""
    Taylor(Ty,x₀,x)

`Ty::∂ℝ` has partials to arbitrary order evaluated at `x₀`. These
partials define a Taylor expansion, which `Taylor` evaluates at value `x`

`Taylor` handles nested structures of `Tuple`s and `SVector`s of `∂ℝ`, applying the
expansion to each element.

See also: [`compose`](@ref), [`McLaurin`](@ref), [`revariate`](@ref), [`fast`](@ref)    
"""
Taylor(y::Tuple,x₀,x) = McLaurin(y,x-x₀)

"""
    compose(Ty,x)

Compose automatic differentiation.  For example 
    Tx = revariate(x)
    Ty = f(Tx)
    y  = compose(Ty,x)    
is faster than
    y  = f(x)
if the length of `x` is smaller than the length of its partials.

See also: [`revariate`](@ref), [`fast`](@ref)    
"""
compose(Ty,x) = McLaurin(Ty,x-VALUE(x))

"""
    y,... = fast(f,x)

In the context of forward automatic differentiation using `∂ℝ`, accelerate the evaluation of
`y,...= f(x)` if the length of `x` is smaller than the length of its partials.

Be extremely careful with closures, making sure that `f` does not capture variables of type `∂ℝ`.

Wrapper function of [`revariate`](@ref) and [`McLaurin`](@ref)      
"""
fast(      f,x) = compose(f(revariate(x)),x)    
justinvoke(f,x) = f(x)    

"""
    composevalue{P,ND}(Ty,X_)

Given `Ty` obtained using `revariate`, and `X_`, obtained using `motion{P}(X)` where `X` is a tuple
of length `ND` and `P=constants(X)`, compute `y`, a tuple of length `ND` of `AbstractArrays` of same `eltype` as vectors in `X.

See also [`revariate`](@ref), [`motion`](@ref), [`motion⁻¹`](@ref), [`composeJacobian`](@ref)  
"""
struct composevalue{P,ND} end
composevalue{P,ND}(Ty,X_) where{P,ND} = motion⁻¹{P,ND}(compose(value{P}(Ty),X_))
composevalue{P,ND}(Ty::Union{Tuple,NamedTuple},X₀) where{P,ND} = map(Tyᵢ->value{P,ND}(Tyᵢ,X₀),Ty)
"""
    composeJacobian{P}(Ty,X_)

Given `Ty` obtained using `revariate`, and `X_`, obtained using `motion{P}(X)` where `X` is a tuple
of `SVectors` and `P=constants(X)`, compute `y`, a tuple of length `ND` of `AbstractArrays` of same `eltype` as vectors in `X,
and `y∂X₀`, the Jacobian of `∂0(y)` with respect to `∂0(X)`.

See also [`revariate`](@ref), [`motion`](@ref), [`motion⁻¹`](@ref), [`composevalue`](@ref)   
"""
struct composeJacobian{P} end
composeJacobian{P}(Ty,X₀) where{P} = compose(∂{P,npartial(Ty)}(Ty),X₀) # y∂X₀
composeJacobian{P}(Ty::Union{Tuple,NamedTuple},X₀) where{P} = map(Tyᵢ->composeJacobian{P}(Tyᵢ,X₀),Ty)

# ∂ℝ( ∂ℝ(a,aₓ), ∂ℝ(aₓ,aₓₓ) ) → ∂ℝ(a,aₓ)   
firstorderonly(a...;)            = firstorderonly.(a)
firstorderonly(a::Tuple)         = firstorderonly.(a)
firstorderonly(a::AbstractArray) = firstorderonly.(a)
firstorderonly(a::∂ℝ)            = precedence(a)≤1 ? a : firstorderonly(a.x) 
firstorderonly(a)                = a

# ∂ℝ(a,aₓ) → ∂ℝ( ∂ℝ(a,aₓ), ∂ℝ(aₓ,0) ) 
backtohigherorder(a::SVector{Na,T},::Type{T}) where{T,Na} = a
backtohigherorder(a::SVector{Na,∂ℝ{1,N,𝕣}},::Type{∂ℝ{2,N,∂ℝ{1,N,𝕣}}}) where{N,Na} = 
     SV{Na}(∂ℝ{2,N,∂ℝ{1,N,𝕣}}( 
                              a[ia],  
                              SV{N}(∂ℝ{1,N,𝕣}(
                                              a[ia].dx[i],
                                              SV{N,𝕣}(zero(𝕣) for j=1:N)
                                              ) for i=1:N)
                             ) for ia=1:Na)
