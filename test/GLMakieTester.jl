struct SpyAxe
    call::Vector{Any}
end
SpyAxe() = SpyAxe(Any[])
lines!(axe,args...;kwargs...)   = push!(axe.call,(fun=:lines!,args=args,kwargs=kwargs))
scatter!(axe,args...;kwargs...) = push!(axe.call,(fun=:scatter!,args=args,kwargs=kwargs))