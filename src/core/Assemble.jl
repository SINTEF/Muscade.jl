# TODO consider Yota.jl

####### Lagrangian from residual and residual from Lagrangian
# an assembler that calls "lagrangian" will call the element's own method if implemented, or this one, which then calls the element's residual method
lagrangian(eleobj::E,δX,X,U,A, t,ε,dbg) where{E<:AbstractElement} = δX ∘₁ residual(eleobj,X,U,A, t,ε,dbg)
# an assembler that calls "residual" will call the element's own method if implemented, or this one, which then calls the element's lagrangian method
function residual(eleobj::E, X,U,A, t,ε,dbg) where{E<:AbstractElement} 
    P            = constants(∂0(X),∂0(U),A,t)
    Nx           = length(∂0(X))
    δX           = δ{P,Nx,𝕣}()                        
    L            = lagrangian(eleobj,δX,X,U,A, t,ε,dbg)
    return ∂{P,Nx}(L)
end
# if an element implements neither lagrangian nor residual, the above code will flat-spin recursively

####### For testing: get all the gradients. 
function gradient(eleobj::E,Λ,X,U,A, t,ε,dbg) where{E<:AbstractElement}
    P            = constants(Λ,∂0(X),∂0(U),A,t)
    nX,nU,nA     = length(Λ),length(∂0(U)),length(A)
    N            = 2nX+nU+nA
    iΛ,iX,iU,iA  = (1:nX) , (1:nX) .+ nX , (1:nU) .+ 2nX , (1:nA) .+ (2nX+nU)  
    ΔY           = δ{P,N,𝕣}()                        
    L            = lagrangian(eleobj,Λ+ΔY[iΛ],(∂0(X)+ΔY[iX],),(∂0(U)+ΔY[iU],),A+ΔY[iA], t,ε,dbg)
    Ly           = ∂{P,N}(L)
    return (L=value{P}(L), Lλ=Ly[iΛ], Lx=Ly[iX], Lu=Ly[iU], La=Ly[iA])
end

###### scaling functions
function scaledlagrangian(scale,eleobj::E,Λs,Xs,Us,As, t,ε,dbg) where{E<:AbstractElement}
    Λ     =       Λs.*scale.Λ                 
    X     = Tuple(xs.*scale.X for xs∈Xs)
    U     = Tuple(us.*scale.U for us∈Us)
    A     =       As.*scale.A
    return lagrangian(eleobj,Λe,Xe,Ue,Ae, t,ε,dbg)
end    
function scaledresidual(scale,eleobj::E, Xs,Us,As, t,ε,dbg) where{E<:AbstractElement} 
    X     = Tuple(xs.*scale.X for xs∈Xs)
    U     = Tuple(us.*scale.U for us∈Us)
    A     =       As.*scale.A
    Re    = scale.Λ .* residual(eleobj, X,U,A, t,ε,dbg) 
end

######## The disassembler
copies(n,a::T) where{T}    = NTuple{n,T}(deepcopy(a) for i∈1:n) # TODO move to Dialect.jl
# dis[ieletyp][iele].index.[  X|U|A]
# dis[ieletyp][iele].scale.[Λ|X|U|A]
struct XUA{T,nX,nU,nA} 
    X::SVector{nX,T}
    U::SVector{nU,T}
    A::SVector{nA,T}
end
struct ΛXUA{T,nX,nU,nA} 
    Λ::SVector{nX,T}
    X::SVector{nX,T}
    U::SVector{nU,T}
    A::SVector{nA,T}
end
struct IS{nX,nU,nA} 
    index:: XUA{𝕫,nX,nU,nA}
    scale::ΛXUA{𝕣,nX,nU,nA}
end
function Disassembler(model::Model)
    neletyp          = length(model.eleobj)  
    dis              = Vector{Any}(undef,neletyp)
    for ieletyp      = 1:neletyp
        nele         = length(model.eleobj[ieletyp])  
        E            = eltype(model.eleobj[ieletyp])
        nX,nU,nA     = getndofs(E)
        dis[ieletyp] = Vector{IS{nX,nU,nA}}(undef,nele)
        iX,iU,iA     =              𝕫1(undef,nX),𝕫1(undef,nU),𝕫1(undef,nA)  # tmp arrays
        sΛ,sX,sU,sA  = 𝕣1(undef,nX),𝕣1(undef,nX),𝕣1(undef,nU),𝕣1(undef,nA)
        for iele     = 1:nele
            ixdof,iudof,iadof = 0,0,0
            for dofID         ∈ model.ele[ieletyp][iele].dofID
                doftyp        = getdoftyp(model,dofID)
                class,scale   = doftyp.class,doftyp.scale
                if     class == :X
                    ixdof    += 1
                    iX[ixdof] = dofID.idof  
                    sX[ixdof] = scale
                    sΛ[ixdof] = scale * model.Λscale
                elseif class == :U
                    iudof    += 1
                    iU[iudof] = dofID.idof
                    sU[iudof] = scale
                elseif class == :A
                    iadof    += 1
                    iA[iadof] = dofID.idof
                    sA[iadof] = scale
                else
                    muscadeerror("dof class must be :X,:U or :A")
                end
            end
            dis[ieletyp][iele] = IS(XUA{𝕫,nX,nU,nA}(iX,iU,iA),ΛXUA{𝕣,nX,nU,nA}(sΛ,sX,sU,sA))
        end
    end
    return dis
end

######## state and initstate
# at each step, contains the complete, unscaled state of the system
struct State{Nxder,Nuder}
    Λ :: 𝕣1
    X :: NTuple{Nxder,𝕣1}
    U :: NTuple{Nuder,𝕣1}
    A :: 𝕣1
    t :: 𝕣
end
# a constructor that provides an initial state
State(model;t=-∞) = State(zeros(getndof(model,:X)),(zeros(getndof(model,:X)),),(zeros(getndof(model,:U)),),zeros(getndof(model,:A)),t)



######## Generic assembler

abstract type Assembler end
function assemble!(asm::Assembler,dis,model,state,ε,dbg)
    zero!(asm)
    for ieletyp ∈ eachindex(model.eleobj)
        eleobj  = model.eleobj[ieletyp]
        assemblesequential!(asm,ieletyp,dis[ieletyp], eleobj,state,ε,(dbg...,ieletyp=ieletyp))
    end
end
function assemblesequential!(asm::Assembler,ieletyp,dis, eleobj,state,ε,dbg) 
    for iele  ∈ eachindex(eleobj)
        scale = dis[iele].scale  # TODO unnecessary replication of "scale": is identical over iele...
        index = dis[iele].index
        Λe    = state.Λ[index.X]                 
        Xe    = Tuple(x[index.X] for x∈state.X)
        Ue    = Tuple(u[index.U] for u∈state.U)
        Ae    = state.A[index.A]
        addin!(asm,scale,ieletyp,iele,eleobj[iele],Λe,Xe,Ue,Ae, state.t,ε,(dbg...,iele=iele))
    end
end

######### generic solver with error management
# TODO move to Dialect.jl
# if a function f is given the argument pointer= Ref{SomeType}()
# the function can then do e.g. vec=allocate(pointer,Vector...) and write to vec.
# and the caller retrievs the data with vec = pointer[] 
# advantage over "return vec" is if f throws, then vec still contains some data.
function allocate(pointer::Ref,target)
    pointer[]=target
    return target
end
function step!(solver!::Function;verbose::𝕓=true,kwargs...) # e.g. solve(SOLstaticX,model,time=1:10)
    verbose && printstyled("\n\n\nMuscade\n",bold=true,color=:cyan)
    pstate = Ref{Any}()
    dbg    = ()
    try
        solver!(pstate,dbg;verbose=verbose,kwargs...) # 
    catch exn
        verbose && report(exn)
        verbose && printstyled("\nAborting the analysis.",color=:red)
        verbose && println(" Function `solve` should still be returning results obtained so far.")
    end
    verbose && printstyled("\nMuscade done.\n\n\n",bold=true,color=:cyan)
    return pstate[]
end

###### DofGroups
abstract type DofGroup end
struct AllXdofs <: DofGroup # TODO add dofgroup reordering here
    scale :: 𝕣1
end
function AllXdofs(model::Model,dis)
    scale  = Vector{𝕣}(undef,getndof(model,:X))
    for di ∈ dis
        for d ∈ di
            scale[d.index.X] = d.scale.X
        end
    end
    return AllXdofs(scale)
end
function Base.setindex!(s::State,x::𝕣1,gr::AllXdofs) # TODO add handling of time derivatives here
    s.X[1] .= x.*gr.scale
end
Base.getindex(s::State,gr::AllXdofs) = s.X[1]./gr.scale
getndof(gr::AllXdofs) = length(gr.scale)

##### Solvers and their Addin
# NB: A solver may require several Assemblers.  Assemblers are object, solvers are functions.

# ASMstaticX: for good old static FEM
struct ASMstaticX <: Assembler 
    dis   :: Vector{Any}          # naïve version! 
    R     :: 𝕣1
    K     :: SparseMatrixCSC{𝕣,𝕫} 
end #  
function ASMstaticX(model::Model,dis) 
    nX       = getndof(model,:X)
    return ASMstaticX(dis,zeros(nX),sparse(Int64[],Int64[],Float64[],nX,nX))
end
function zero!(asm::ASMstaticX)
    asm.R  .= 0
    asm.K  .= 0
end
function addin!(asm::ASMstaticX,scale,ieletyp,iele,eleobj::E,Λ,X,U,A, t,ε,dbg)  where{E<:AbstractElement}
    Nx           = length(Λ) 
    ΔX           = δ{1,Nx,𝕣}()                 # NB: precedence==1, input must not be Adiff
    Re           = scaledresidual(scale,eleobj, (∂0(X)+ΔX,),U,A, t,ε,dbg)
    i            = Vector(asm.dis[ieletyp][iele].index.X)    # TODO not type stable (X is SVector).  Allocating!
    asm.R[i  ]  += value{1}(Re)            
    asm.K[i,i]  += ∂{1,Nx}(Re)                     # TODO very slow!   TODO can a sparse be indexed by a view? or do I need a i-buffer in asm?
end
function StaticX(pstate,dbg;model::Model,time::AbstractVector{𝕣},
                    initial::State=State(model),
                    maxiter::ℤ=50,maxΔy::ℝ=1e-5,maxR::ℝ=∞,
                    verbose::𝕓=true,saveiter::𝔹=false)
    # important: this code assumes that there is no χ in state.
    verb             = verbose
    verb && @printf "    StaticX solver\n\n"
    dis              = Disassembler(model)
    asm              = ASMstaticX(model,dis)
    dofgr            = AllXdofs(model,dis)
    asmt,solt,citer  = 0.,0.,0
    cΔy²,cR²         = maxΔy^2,maxR^2
    state            = allocate(pstate,Vector{State}(undef,saveiter ? maxiter : length(time))) # state is not a return argument so that data is not lost in case of exception
    for (it,t)       ∈ enumerate(time)
        verb && @printf "    increment %3d" it
        old          = it==1 ? initial : state[it-1]
        s            = State(old.Λ,old.X,old.U,old.A,t)
        y            = s[dofgr] # includes scaling
        for iiter    = 1:maxiter
            citer   += 1
            asmt+=@elapsed assemble!(asm,dis,model,s, 0.,(dbg...,solver=:StaticX,it=it,iiter=iiter))
            solt+=@elapsed Δy = try asm.K\-asm.R catch; muscadeerror(@sprintf("Incremental solution failed at it=%i, iiter=%i",it,iiter)) end
            Δy²,R²   = sum(Δy.^2),sum(asm.R.^2)
            y      .+= Δy
            s[dofgr] = y  # includes descaling
            saveiter && (state[iiter]=s)
            if Δy²≤cΔy² && R²≤cR² 
                verb && @printf " converged in %3d iterations. |Δy|=%7.1e |R|=%7.1e\n" iiter √(Δy²) √(R²)
                saveiter || (state[it]=s)
                break#out of the iiter loop
            end
            iiter==maxiter && muscadeerror(@sprintf(" no convergence after %3d iterations |Δy|:%g / %g, |R|:%g / %g",iiter,√(Δy²),maxΔy,√(R²)^2,maxR))
        end
    end
    verb && @printf "\n    nel=%d, ndof=%d, nincr=%d, niter=%d, niter/nincr=%5.2f\n" getnele(model) getndof(dofgr) length(time) citer citer/length(time)
    verb && @printf "    Build  time = %s, (per iteration: %s, per iteration and element: %s)\n" showtime(asmt)  showtime(asmt/citer)  showtime(asmt/citer/getnele(model))
    verb && @printf "    Solve  time = %s, (per iteration: %s, per iteration and dof:     %s)\n" showtime(solt)  showtime(solt/citer)  showtime(solt/citer/getndof(dofgr))
end
