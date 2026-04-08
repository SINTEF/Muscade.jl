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
  by `motion⁻¹`. Touched variables can for example be marked with an underscore.
- The precendence `P` must be calculated using `constants` with all variables that are input to 
  the function and may be differentiated.
- If other levels of automatic differentiation are introduced within the function, unpack in reverse
  order of packing.    

See [`motion⁻¹`](@ref)
"""
motion{ P    }(a::NTuple{ND,SV{N,R}}) where{ND,P,N,R        } = SV{N}(motion_{P-1,P+ND-2}(let a=a 
                                                                                            ntuple(j->a[j][i],Val(ND))
                                                                                         end) for i=1:N)
#motion{ P    }(a::NTuple{ND,SV{N,R}}) where{ND,P,N,R        } = SV{N}(motion_{P-1,P+ND-2}(ntuple(j->a[j][i],Val(ND))) for i=1:N)
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
motion⁻¹{P,1,0}(a::ℝ) where{P} =                      a
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
motion⁻¹{P,ND   }(a::Tuple          ) where{P,ND} = (motion⁻¹{P,ND}(first(a)),motion⁻¹{P,ND}(Base.tail(a))...)   
motion⁻¹{P,ND   }(a::Tuple{}        ) where{P,ND} = ()   
motion⁻¹{P,ND   }(a::NamedTuple     ) where{P,ND} = NamedTuple{keys(a)}(motion⁻¹{P,ND}(values(a)))  
motion⁻¹{P,ND   }(a...              ) where{P,ND} = motion⁻¹{P,ND}(a) 

#############

flat_length(a::NamedTuple   )                = flat_length(values(a))
flat_length(a::Tuple        )                = flat_length(first(a))+flat_length(Base.tail(a))
flat_length(a::Tuple{}      )                = 0
flat_length(a::SArray       )                = length(a) 
flat_length(a::ℝ            )                = 1

flat_eltype(a::NamedTuple   )                = flat_eltype(values(a))
flat_eltype(a::Tuple        )                = promote_type(flat_eltype(first(a)),flat_eltype(Base.tail(a)))
flat_eltype(a::Tuple{X}     ) where{X}       = flat_eltype(a[1]) 
flat_eltype(a::SArray{S,T}  ) where{S,T}     = T  
flat_eltype(a::R            ) where{R<:ℝ}    = R

precedence(a::NamedTuple   )                 = precedence(values(a))
precedence(a::Tuple        )                 = max(precedence(first(a)),precedence(Base.tail(a)))
precedence(a::Tuple{X}     ) where{X}        = precedence(a[1]) 

npartial(  a::NamedTuple   )                 = npartial(values(a))
npartial(  a::Tuple        )                 = max(npartial(first(a)),npartial(Base.tail(a)))
npartial(  a::Tuple{X}     ) where{X}        = npartial(a[1]) 

struct flatten{T} end
flatten(a)                                   = flatten{flat_eltype(a)}(a)
flatten{T}(a::NamedTuple   ) where{T}        = flatten{T}(values(a))
flatten{T}(a::Tuple        ) where{T}        = SVector{flat_length(a),T}(flatten{T}(first(a))..., flatten{T}(Base.tail(a))...)  
flatten{T}(a::Tuple{}      ) where{T}        = ()
flatten{T}(a::SArray       ) where{T}        = SVector{length(a),T}(a)  
flatten{T}(a::ℝ            ) where{T}        = T(a)

struct type_multivariate_𝕣{P,N} 
    dummy::𝕫  
end
type_multivariate_𝕣{0,N}() where{  N} = 𝕣
type_multivariate_𝕣{P,N}() where{P,N} = ∂ℝ{P,N,type_multivariate_𝕣{P-1,N}()} # this causes (slight) type instability - because if Ra isnot ℝ, then return type is different.

struct multivariate_𝕣{P,N} end
multivariate_𝕣{P,N}(a   ,i      ) where{P,N}          = ∂ℝ{P,N  }(multivariate_𝕣{P-1,N}(a,i      ),i) 
multivariate_𝕣{0,N}(a::𝕣,i      ) where{  N}          = a

multivariate_𝕣{P,N}(a::R,i,scale) where{P,N,R}        = ∂ℝ{P,N  }(multivariate_𝕣{P-1,N}(a,i,scale),SV{N,R}(i==j ? scale : zero(R) for j=1:N)) 
multivariate_𝕣{0,N}(a::𝕣,i,scale) where{  N}          = a

""" 
    TX = revariate{P}(V)

The variable `V` is a nested structure `NamedTuple`s, `Tuple`s and `SArrays` of 
`Real`s (possibly: `∂ℝ`s).

`V` is stripped of its partials (if any), an [re]variated to 
order `P`, with npartial `N` equal the "flat length" of `V`.   

    TV = revariate(VX)

revariates to the order `precedence(V)`.  

Be extremely careful never to mix any variable that is a function of `V` with any other variables
containing  `∂ℝ`s but not produced by the same `revariate`.

# Optimizing performance of automatic differentiation

`revariate`, in conjunction with `chainrule` can be used to improve performance when the length of 
`V` is smaller than the length of its partials.

See also: [`chainrule`](@ref), [`Taylor`](@ref), [`McLaurin`](@ref), [`fast`](@ref) 

# Relation to `variate`

`revariate` differs from [`variate`](@ref) in two ways:
1) `variate` does not strip its input of partials. It wraps the input with new partials.
2) `variate` can only handle `SVector` inputs.  `revariate` can handle a broader variety of inputs
   which makes it convenient to differentiate with resepct to a collection of `SVector`.

# Scaling

    V = (;X,U,A)
    S = (X=scale.X,U=scale.U,A=scale.A)
    TV = revariate{P}(V,S)

allows to introduce scaling in automatic differentiation.  For this method, `S` and `V`
have the same structure, with the important exception that `Tuple`s in `V` must be
`ntuple`s (have elements of identical type `T`), and, to a `Tuple` in `V` corresponds 
a variable of type `T` in `V`. Put simply: the same scale `scale.X` will be applied
to `∂0[X]`, `∂1[X]` and `∂2[X]`. 

See also: [`chainrule`](@ref), [`to_order`](@ref)
"""
struct revariate{P,N,Z}   end
revariate(a)                                               = revariate{precedence(a)}(a)
revariate{P           }(a                 ) where{P}       = revariate{P,flat_length(a),:variate}(a,1)
revariate{P,N,Z       }(a::NamedTuple   ,i) where{P,N,Z}   = NamedTuple{keys(a)}(revariate{P,N,Z}(values(a),i)) 
revariate{P,N,Z       }(a::Tuple        ,i) where{P,N,Z}   = (revariate{P,N,Z}(first(a),i),revariate{P,N,Z}(Base.tail(a),i+flat_length(first(a)))...)
revariate{P,N,Z       }(a::Tuple{}      ,i) where{P,N,Z}   = ()
revariate{P,N,Z       }(a::SArray{S}    ,i) where{P,N,Z,S} = SArray{S,type_multivariate_𝕣{P,N}()}(revariate{P,N,Z}(aⱼ,i-1+j) for (j,aⱼ)∈enumerate(a))
revariate{P,N,:δ      }(a::ℝ            ,i) where{P,N  }   = multivariate_𝕣{P,N}(zero( a),i)
revariate{P,N,:variate}(a::ℝ            ,i) where{P,N  }   = multivariate_𝕣{P,N}(VALUE(a),i)

# Specialised to adiff dofs, with scaling
# s for scale
# IMPORTANT: in a Tuple (but not NamedTuple), all a[i] share the same s
struct revariatevalues{P,N,Z}   end
revariate(a,s)                                                            = revariate{precedence(a)}(a,s)
revariate{P           }(a                 ,s             ) where{P}       = revariate{P,flat_length(a),:variate}(a,1,s)
revariate{P,N,Z       }(a::NamedTuple   ,i,s             ) where{P,N,Z}   = NamedTuple{keys(a)}(revariatevalues{P,N,Z}(values(a),i,values(s))) 
revariatevalues{P,N,Z }(a::Tuple        ,i,s             ) where{P,N,Z}   = (revariate{P,N,Z}(first(a),i,first(s)),revariatevalues{P,N,Z}(Base.tail(a),i+flat_length(first(a)),Base.tail(s))...)
revariatevalues{P,N,Z }(a::Tuple{}      ,i,s             ) where{P,N,Z}   = ()
revariate{P,N,Z       }(a::Tuple        ,i,s             ) where{P,N,Z}   = (revariate{P,N,Z}(first(a),i,      s ),revariate{      P,N,Z}(Base.tail(a),i+flat_length(first(a)),          s )...)
revariate{P,N,Z       }(a::Tuple{}      ,i,s             ) where{P,N,Z}   = ()
revariate{P,N,Z       }(a::SArray{S}    ,i,s::SArray{S,𝕣}) where{P,N,Z,S} = SArray{S,type_multivariate_𝕣{P,N}()}(revariate{P,N,Z}(a[j],i-1+j,s[j]) for j∈eachindex(a))
revariate{P,N,:δ      }(a::ℝ            ,i,s::𝕣          ) where{P,N  }   = multivariate_𝕣{P,N}(zero( a),i,s)
revariate{P,N,:variate}(a::ℝ            ,i,s::𝕣          ) where{P,N  }   = multivariate_𝕣{P,N}(VALUE(a),i,s)

# create an "increment" δ (same as revariate, but value set to zero)
struct reδ{P} end
reδ{P}(a  ) where{P} = revariate{P,flat_length(a),:δ}(a,1  )
reδ{P}(a,s) where{P} = revariate{P,flat_length(a),:δ}(a,1,s)

# works with revariate{P}(a,s) and returns a structure of indices into the partials
revariate_indices( args...)      = (revariate_indices_(args,0)...,flat_length(args))
revariate_indices_(a::Tuple  ,i) = (revariate_indices_(first(a),i),revariate_indices_(Base.tail(a),i+flat_length(first(a)))...)
revariate_indices_(a::Tuple{},i) = ()
revariate_indices_(a::SArray ,i) = i+1:i+flat_length(a)
revariate_indices_(a::𝕣,i)       = i+1

"""
    to_order{P}(V)

Decrease (lossy) or increase (pad partials with zeros) the order of differentiation of `V`.
`V` is a nested structure of `NamedTuple`, `Tuple`, `SArray`, and the components
of `V` must be of type `∂ℝ` (otherwise, the number of partials would be undefined).
"""
struct to_order{P,N} end
to_order{P,N}(a::NamedTuple   ) where{P,N} = NamedTuple{keys(a)}(to_order{P,N}(values(a)))
to_order{P,N}(a::Tuple        ) where{P,N} = (to_order{P,N}(first(a)),to_order{P,N}(Base.tail(a))...)
to_order{P,N}(a::Tuple{}      ) where{P,N} = ()
to_order{P,N}(a::AbstractArray) where{P,N} = to_order{P,N}.(a)
function to_order{P,N}(a::Ra)   where{Ra<:∂ℝ{Pa,N},P} where{Pa,N}
    if     Pa==P
        a
    elseif Pa> P   
        to_order{P,N}(value{Pa}(a))
    elseif Pa< P
        ap  = to_order{P-1,N}(a) 
        ∂ℝ{P,N}(ap, SVector{N}(∂ℝ{P-1,N}(ap.dx[i]) for i=1:N) )
    end
end
function to_order{P,N}(a::𝕣) where{P,N} 
    P==0 ? a : ∂ℝ{P,N}(to_order{P-1,N}(a))
end
to_order{Pa}(a) where{Pa}= to_order{Pa,npartial(a)}(a)
 



#########################

"""
    McLaurin(Ty,x)

`Ty::∂ℝ` has partials to arbitrary order with respect to a variable `x`. These
partials define a McLaurin expansion, which `McLaurin` evaluates at value `x`, 
as if `Ty` had been computed at 0.

`McLaurin` handles nested structures of `Tuple`s and `SVector`s of `∂ℝ`, applying the
expansion to each element.

`McLaurin` is a utility function behind [`chainrule`](@ref) and [`Taylor`](@ref)

See also: [`chainrule`](@ref), [`Taylor`](@ref), [`revariate`](@ref), [`fast`](@ref)    
"""
McLaurin(y::Tuple,Δx)                          = tuple(McLaurin(first(y),Δx),McLaurin(Base.tail(y),Δx)...) 
McLaurin( ::Tuple{},Δx)                        = tuple() 
McLaurin(y::SArray{Sy,Ty,Dy,Ly},Δx::SVector{Sx,Tx}) where{Sy,Ty,Dy,Ly,Sx,Tx} = SArray{Sy,Tx,Dy,Ly}(McLaurin(yᵢ,Δx) for yᵢ∈y)
McLaurin(y::SArray{Sy,𝕣 ,Dy,Ly},Δx::SVector{Sx,Tx}) where{Sy,   Dy,Ly,Sx,Tx} =                              y
McLaurin(y::∂ℝ,Δx)                             = McLaurin(y.x,Δx) .+ McLaurin_right(y,Δx)
McLaurin(y::𝕣 ,Δx)                             =          y
function McLaurin_right(y::∂ℝ{P,N,R},Δx::SVector{N}) where{P,N,R} 
    if N==0
        return zero(y) # hum!!!!
    else
        s = McLaurin_right(y.dx[1],Δx)*Δx[1]
        for i ∈ 2:N
            s += McLaurin_right(y.dx[i],Δx)*Δx[i]
        end
        s /= P
    end
    return s
end
McLaurin_right(y::𝕣    ,Δx            )        =          y

"""
    Taylor(Ty,x₀,x)

`Ty::∂ℝ` has partials to arbitrary order evaluated at `x₀`. These
partials define a Taylor expansion, which `Taylor` evaluates at value `x`

`Taylor` handles nested structures of `Tuple`s and `SVector`s of `∂ℝ`, applying the
expansion to each element.

See also: [`chainrule`](@ref), [`McLaurin`](@ref), [`revariate`](@ref), [`fast`](@ref)    
"""
Taylor(y::Tuple,x₀,x) = McLaurin(y,x-x₀)

"""
    chainrule(Ty,x)

Apply a chain rule in automatic differentiation.  For example 
    Tx = revariate(x)
    Ty = f(Tx)
    y  = chainrule(Ty,x)    
is faster than
    y  = f(x)
if the length of `x` is smaller than the length of its partials.

See also: [`revariate`](@ref), [`fast`](@ref)    
"""
function chainrule(Ty,x) 
    fx = flatten(x)
    McLaurin(Ty,fx.-VALUE(fx))
end

"""
    y,... = fast(f,x)

In the context of forward automatic differentiation using `∂ℝ`, accelerate the evaluation of
`y,...= f(x)` if the length of `x` is smaller than the length of its partials.

Also work where `x` is a nested structure of `Tuple`s and `NamedTuple`s where the leaves
are `ℝ` or `SArray{S,R} where {S,R<:ℝ}`.    

Be extremely careful with closures, making sure that `f` does not capture variables of type `∂ℝ`.

Wrapper function of [`revariate`](@ref) and [`McLaurin`](@ref)      
"""
fast(      f,x) = apply{:chainrule}(f,x)    
justinvoke(f,x) = apply{:direct   }(f,x)    


struct apply{Mode} end
apply{:chainrule}(f,x) = chainrule(f(revariate(x)),x)
apply{:direct   }(f,x) = f(x)

"""
    composevalue{P,ND}(Ty,X_)

Given `Ty` obtained using `revariate`, and `X_`, obtained using `motion{P}(X)` where `X` is a tuple
of length `ND` and `P=constants(X)`, compute `y`, a tuple of length `ND` of `AbstractArrays` of same `eltype` as vectors in `X.

See also [`revariate`](@ref), [`motion`](@ref), [`motion⁻¹`](@ref), [`composeJacobian`](@ref)  
"""
struct composevalue{P,ND} end
composevalue{P,ND}(Ty            ,X_) where{P,ND} = motion⁻¹{P,ND}(chainrule(value{P}(Ty),X_))
composevalue{P,ND}(Ty::NamedTuple,X_) where{P,ND} = NamedTuple{keys(Ty)}(composevalue{P,ND}(values(Ty),X_))
composevalue{P,ND}(Ty::Tuple     ,X_) where{P,ND} = (composevalue{P,ND}(first(Ty),X_),composevalue{P,ND}(Base.tail(Ty),X_)...)
composevalue{P,ND}(Ty::Tuple{}   ,X_) where{P,ND} = ()
"""
    composeJacobian{P}(Ty,X_)

Given `Ty` obtained using `revariate`, and `X_`, obtained using `motion{P}(X)` where `X` is a tuple
of `SVectors` and `P=constants(X)`, compute `y`, a tuple of length `ND` of `AbstractArrays` of same `eltype` as vectors in `X,
and `y∂X₀`, the Jacobian of `∂0(y)` with respect to `∂0(X)`.

See also [`revariate`](@ref), [`motion`](@ref), [`motion⁻¹`](@ref), [`composevalue`](@ref)   
"""
struct composeJacobian{P} end
composeJacobian{P}(Ty            ,X₀) where{P} = chainrule(∂{P,npartial(Ty)}(Ty),X₀) # y∂X₀
composeJacobian{P}(Ty::NamedTuple,X₀) where{P} = NamedTuple{keys(Ty)}(composeJacobian{P}(values(Ty),X₀))
composeJacobian{P}(Ty::Tuple     ,X₀) where{P} = (composeJacobian{P}(first(Ty),X₀),composeJacobian{P}(Base.tail(Ty),X₀)...)
composeJacobian{P}(Ty::Tuple{}   ,X₀) where{P} = ()


