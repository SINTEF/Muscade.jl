struct SpyAxe
    call::Vector{Any}
end
SpyAxe() = SpyAxe(Any[])
lines!(axe::SpyAxe,args...;kwargs...)   = push!(axe.call,(fun=:lines!,args=args,kwargs=kwargs))
scatter!(axe::SpyAxe,args...;kwargs...) = push!(axe.call,(fun=:scatter!,args=args,kwargs=kwargs))
mesh!(axe::SpyAxe,args...;kwargs...) = push!(axe.call,(fun=:mesh!,args=args,kwargs=kwargs))

