struct SpyAxe
    data::Vector{Any}
end
SpyAxe() = SpyAxe(Any[])
lines!(axe,args...;kwargs...)   = push!(axe.data,(fun=:lines!,args=args,kwargs=kwargs))
scatter!(axe,args...;kwargs...) = push!(axe.data,(fun=:scatter!,args=args,kwargs=kwargs))