
struct Taylor{O,Nx,TA}
    x::SVector{Nx,𝕣}
    A::TA
end
Taylor{O}(x::SVector{Nx,𝕣},A::TA) where{O,Nx,TA} = Taylor{O,Nx,TA}(x,A)
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
Taylor(f::Function,X::SVector{Nx,R}) where{Nx,R<:Real} = Taylor{min(precedence(R),2)}(f,X) 
function Taylor{0}(f::Function,X::SVector{Nx,R}) where{Nx,R<:Real}
    x  = VALUE(X)
    y  = f(x)
    A  = (y,)
    return Taylor{0}(x,A)
end    
function Taylor{1}(f::Function,X::SVector{Nx,R}) where{Nx,R<:Real}
    x  = VALUE(X)
    y  = f(variate{1,Nx}(x))
    A  = value_∂{1,Nx}(y)  # 
    return Taylor{1}(x,A)
end    
function Taylor{2}(f::Function,X::SVector{Nx,R}) where{Nx,R<:Real}
    x  = VALUE(X)
    y  = f(variate{2,Nx}(variate{1,Nx}(x)))
    A  = (value{1}(value{2}(y)), 
          ∂{1,Nx}( value{2}(y)), 
          ∂{1,Nx}(∂{2,Nx}(y))/2)
    return Taylor{2}(x,A)
end    

(te::Taylor{0,Nx,TA})(X::SVector{Nx,R}) where{Nx,R<:Real,TA} = te.A[1]
(te::Taylor{1,Nx,TA})(X::SVector{Nx,R}) where{Nx,R<:Real,TA} = (te.A[2])∘₁(X-te.x).+te.A[1]
function (te::Taylor{2,Nx,TA})(X::SVector{Nx,R}) where{Nx,R<:Real,TA}
    dX = X-te.x
    return (te.A[3]∘₁dX .+ te.A[2])∘₁dX.+te.A[1]
end

∂(t::Taylor{O}) where{O} = Taylor{O-1}(t.x,t.A[2:end])
∂(t::Taylor{0}) where{O} = MuscadeError("Tried to differentiate a 0th order Taylor")
