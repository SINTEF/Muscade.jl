
# TODO: FunctionWrappers.jl

"""
    f  = FunctionFromVector(xs::AbstractRange,ys::AbstractVector)
    y  = f(x)

    Linear interpolation.  Fails if `x` is outside the range `xs`
"""
struct FunctionFromVector{X,Y} <:Function
    x::X  
    y::Y
end
FunctionFromVector(x::AbstractRange,y::AbstractVector{R}) where{R<:Real} = FunctionFromVector{typeof(x),typeof(y)}(x,y)
function (f::FunctionFromVector)(x)
    ix = (x-first(f.x))/step(f.x) +1
    @assert 1≤ix≤length(f.x)

    i  = min(trunc(𝕫,ix),length(f.x)-1)
    di = ix-i
    return f.y[i]*(1-di)+f.y[i+1]*di
end

struct QuadraticFunctionWithConstantMean <: Function
    μ::𝕣
    σ::𝕣
end
struct QuadraticFunctionWithMeanFuncOfTime{F} <:Function
    μ::F
    σ::𝕣
end
"""
    f  = QuadraticFunction(μ,σ)

`μ` and `σ` are `𝕣` (`Float64`)

    y  = f(x) # == 1/2*((x-μ)/σ)^2
    f  = QuadraticFunction(μ,σ)

Alternatively, `μ` can be a `Function` of time, in which case

    y  = f(x,t) # == 1/2*((x-μ(t))/σ)^2

"""
QuadraticFunction(μ::𝕣       ,σ::𝕣) = QuadraticFunctionWithConstantMean(  μ,σ)
QuadraticFunction(μ::Function,σ::𝕣) = QuadraticFunctionWithMeanFuncOfTime(μ,σ)

(f::QuadraticFunctionWithConstantMean  )(x,args...) = .5*((x-f.μ   )/f.σ)^2  # args... allows to ignore an extra t argument
(f::QuadraticFunctionWithMeanFuncOfTime)(x,t) = .5*((x-f.μ(t))/f.σ)^2