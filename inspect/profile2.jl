using Profile,ProfileView,BenchmarkTools
using Muscade
using StaticArrays,GLMakie,Random,Printf#nλfloor

#------------- custom element code
struct Pipe <: AbstractElement
    K      :: SMatrix{4,4,𝕣,16}
    W      :: SVector{4,𝕣}
    x      :: SVector{2,𝕣}
end
function Pipe(nod::Vector{Node};EI,w) 
    c = coord(nod)
    x1,x2 = c[1][1],c[2][1]
    L = abs(x2-x1)
    L² = L^2  
    K = EI/L^3*SMatrix{4,4,𝕣}(12,6L,-12,6L, 6L,4L²,-6L,2L², -12,-6L,12,-6L, 6L,2L²,-6L,4L²)
    W = w*L*SVector{4,𝕣}(1,0,1,0)
    return Pipe(K,W,SVector(x1,x2))
end
Muscade.doflist(::Type{Pipe}) =(inod=(1,1,2,2), class=(:X,:X,:X,:X), field=(:z,:r,:z,:r))
@espy function Muscade.residual(o::Pipe, X,U,A, t,γ,dbg) 
    return o.K*∂0(X)+o.W
end
function Muscade.draw(axe,key,out, o::Pipe, δX,X,U,A, t,γ,dbg;kwargs...)
    z    = ∂0(X)[SVector(1,3)]  
    lines!(axe,collect(o.x),collect(z) ;kwargs...)
end
#---------------- analysis and result extraction
nel  = 2000
δx   = 1.
δz   = 0.1
w    = 0#100*9.81
EI   = 5000000
Δz   = 0.25 # 0.25 causes a rattle if not interior point

x               = δx*(0:nel)
z               = zeros(nel+1)
rng             = MersenneTwister(0)
zfloor          = cumsum(0.1*randn(rng,nel+1))
zfloor        .-= range(0,zfloor[end],nel+1) .- .01
@once floorgap(z,t,zfloor)=z[1]-zfloor
@once ceilgap(z,t,zceil)=zceil-z[1]
# zceil           = zfloor.+Δz
# floorgap        = [(z,t)->z[1]-zfloor[i] for i=1:nel+1] # comment to avoid recompilation
# ceilgap         = [(z,t)->zceil[i]-z[1]  for i=1:nel+1] # comment to avoid recompilation
coords          = cat(x,z,dims=2)
gₛ               = 0.1
λₛ               = 4000.   

model           = Model(:Pipe)
node            = addnode!(model,coords) 
pipe            = [addelement!(model,Pipe,node[i:i+1];EI,w) for i=1:nel]
floorctct       = [addelement!(model,Constraint,[node[i]];λₛ,gₛ,xinod=(1,),xfield=(:z,),
                               λinod=1, λclass=:X, λfield=:λfloor, g=floorgap,gargs=(zfloor[i],), mode=inequal) for i=1:nel+1]
ceilctct        = [addelement!(model,Constraint,[node[i]];λₛ,gₛ,xinod=(1,),xfield=(:z,),
                               λinod=1, λclass=:X, λfield=:λceil,  g=ceilgap,gargs=(zfloor[i]+Δz,), mode=inequal) for i=1:nel+1]

initialstate    = initialize!(model)
mission = :profile
if mission == :report
    state           = solve(StaticX;initialstate,time=[0.],maxiter=200,saveiter=true,γ0=10.,γfac1=.5,γfac2=10.,verbose=true)
elseif mission == :time
    state           = solve(StaticX;initialstate,time=[0.],maxiter=200,saveiter=true,γ0=10.,γfac1=.5,γfac2=10.,verbose=false)
    @btime state    = solve(StaticX;initialstate,time=[0.],maxiter=200,saveiter=true,γ0=10.,γfac1=.5,γfac2=10.,verbose=false)
elseif mission == :profile
    Profile.clear()
    Profile.@profile for i=1:10
        local state           = solve(StaticX;initialstate,time=[0.],maxiter=200,saveiter=true,γ0=10.,γfac1=.5,γfac2=10.,verbose=false)
    end
ProfileView.view(fontsize=30);
end


if false # sandwich plot
    state           = cat(dims=1,initialstate,state[1:findlastassigned(state)])
    niter           = length(state)

    out,key         = getresult(state,@request((λ,g)),floorctct)
    λfloor          = out[key.λ,:,:] # out[ikey,iele,istep]
    gfloor          = out[key.g,:,:] # out[ikey,iele,istep]

    out,key         = getresult(state,@request((λ,g)),ceilctct)
    λceil           = out[key.λ,:,:] # out[ikey,iele,istep]
    gceil           = out[key.g,:,:] # out[ikey,iele,istep]

    z,dofid         = getdof(state,field=:z) #x[iel,idir,iiter]

    fig1      = Figure(resolution = (1000,250))
    with_theme(Theme(fontsize = 30,font="Arial")) do
        axe      = Axis(fig1[1,1])
        hidespines!(axe)
        hidedecorations!(axe)
        for i=1:niter-1
            lines!(axe,x,z[:,1,i];linewidth=1,color=:grey)
        end
        draw(axe,state[niter];linewidth=2,color=:black)
        lines!(axe,x,zfloor,color=:red, linewidth=2)
        lines!(axe,x,zceil,color=:red, linewidth=2)
    end
    display(fig1)
end
;