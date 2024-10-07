using GLMakie
struct SpyAxe
    call::Vector{Any}
end
SpyAxe() = SpyAxe(Any[])
GLMakie.lines!(axe::SpyAxe,args...;kwargs...)   = push!(axe.call,(fun=:lines!,args=args,kwargs=kwargs))
GLMakie.scatter!(axe::SpyAxe,args...;kwargs...) = push!(axe.call,(fun=:scatter!,args=args,kwargs=kwargs))

