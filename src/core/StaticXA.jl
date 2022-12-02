###### DofGroups
struct AllAdofs <: DofGroup 
    scale :: 𝕣1
end
function AllAdofs(model::Model,dis)
    scale  = Vector{𝕣}(undef,getndof(model,:A))
    for di ∈ dis
        for d ∈ di
            scale[d.index.A] = d.scale.A
        end
    end
    return AllAdofs(scale)
end
function Base.setindex!(s::State,a::𝕣1,gr::AllAdofs) 
    s.A .= a.*gr.scale
end
Base.getindex(s::State,gr::AllAdofs) = s.A./gr.scale
getndof(gr::AllAdofs) = length(gr.scale)


struct AllΛXUdofs <: DofGroup 
    Λscale :: 𝕣1
    Xscale :: 𝕣1
    Uscale :: 𝕣1
    nX     :: 𝕣
    nU     :: 𝕣
end
function AllΛXUdofs(model::Model,dis)
    nX     = getndof(model,:X)
    nU     = getndof(model,:U)
    Λscale = Vector{𝕣}(undef,nX)
    Xscale = Vector{𝕣}(undef,nX)
    Uscale = Vector{𝕣}(undef,nU)
    for di ∈ dis
        for d ∈ di
            Λscale[d.index.Λ] = d.scale.Λ
            Xscale[d.index.X] = d.scale.X
            Uscale[d.index.U] = d.scale.U
        end
    end
    return AllΛXUdofs(Λscale,Xscale,Uscale,nX,nU)
end
function Base.setindex!(s::State,y::𝕣1,gr::AllΛXUdofs) 
    s.Λ .= y[      1:      s.nX].*gr.Λscale
    s.X .= y[ s.nX+1:     2s.nX].*gr.Xscale
    s.U .= y[2s.nX+1:2s.nX+s.nU].*gr.Uscale
end
function Base.getindex(s::State,gr::AllΛXUdofs) 
    y = 𝕣1(undef,2s.nX+s.nU)
    y[      1:      s.nX] = s.Λ ./gr.Λscale
    y[ s.nX+1:     2s.nX] = s.X ./gr.Xscale
    y[2s.nX+1:2s.nX+s.nU] = s.U ./gr.Uscale
    return y
end
getndof(gr::AllΛXUdofs) = length(2gr.nX+gr.nU)


##### Solvers and their Addin
# NB: A solver may require several Assemblers.  Assemblers are object, solvers are functions.

# ASMstaticX: for good old static FEM
struct ASMstaticΛXU_A <: Assembler 
    dis   :: Vector{Any}          # naïve version! 
    Ly    :: 𝕣1
    La    :: 𝕣1
    Lyy   :: SparseMatrixCSC{𝕣,𝕫} 
    Lya   :: SparseMatrixCSC{𝕣,𝕫} 
    Laa   :: SparseMatrixCSC{𝕣,𝕫} 
    nX    :: 𝕫
end #  
spa(a,b) = sparse(Int64[],Int64[],Float64[],a,n)
function ASMstaticΛXU_A(model::Model,dis) 
    nX     = getndof(model,:X)
    nU     = getndof(model,:U)
    nA     = getndof(model,:A)
    return ASMstaticΛXU_A(dis,zeros(2nX+nU),zeros(nA),spa(2nX+nU,2nX+nU),spa(2nX+nU,nA),spa(nA,nA),nX)
end
function zero!(asm::ASMstaticΛXU_A)
    asm.Ly  .= 0
    asm.La  .= 0
    asm.Lyy .= 0
    asm.Lya .= 0
    asm.Laa .= 0
end
function addin!(asm::ASMstaticΛXU_A,scale,ieletyp,iele,eleobj::E,Λ,X,U,A, t,ε,dbg)  where{E<:AbstractElement}
    Nx           = length(Λ) # in the element
    Nu           = length(U[1])
    Na           = length(A)   
    NX           = asm.nX # in the model                
    Δz           = δ{      1,2Nx+Nu+Na,𝕣}(  )                   # NB: precedence==1, input must not be Adiff 
    ΔZ           = variate{2,2Nx+Nu+Na  }(Δz)
    ΔΛ           = view(ΔZ,       1: Nx      ) # TODO Static?
    ΔX           = view(ΔZ, Nx   +1:2Nx      ) # TODO Static?
    ΔU           = view(ΔZ,2Nx   +1:2Nx+Nu   ) # TODO Static?
    ΔA           = view(ΔZ,2Nx+Nu+1:2Nx+Nu+Na) # TODO Static?
    L            = scaledlagrangian(scale,eleobj, Λ+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA, t,ε,dbg)
    iX           = Vector(asm.dis[ieletyp][iele].index.X)  # sparse doesn't like static indices  
    iU           = Vector(asm.dis[ieletyp][iele].index.U)    
    iA           = Vector(asm.dis[ieletyp][iele].index.A) 
    ∇L           = ∂{2,Nx}(L)
    Lz           = value{1}(∇L)
    Lzz          = ∂{1,Nx}(∇L)
    asm.Ly[iX      ]  += Lz[       1: Nx      ]  # Lλ
    asm.Ly[iX+ NX  ]  += Lz[ Nx   +1:2Nx      ]  # Lx 
    asm.Ly[iU+2NX  ]  += Lz[2Nx   +1:2Nx+Nu   ]  # Lu 
    asm.La[iA      ]  += Lz[2Nx+Nu+1:2Nx+Nu+Na]  # Lu 
              
    asm.K[i,i]  += ∂{1,Nx}(L)                     
end
function StaticX(pstate,dbg;model::Model,time::AbstractVector{𝕣},
                    initial::State=State(model,Disassembler(model)),
                    maxiter::ℤ=50,maxΔy::ℝ=1e-5,maxR::ℝ=∞,
                    verbose::𝕓=true,saveiter::𝔹=false)
    # important: this code assumes that there is no χ in state.
    verb             = verbose
    verb && @printf "    StaticX solver\n\n"
    dis              = initial.dis
    asm              = ASMstaticX(model,dis)
    dofgr            = AllXdofs(model,dis)
    asmt,solt,citer  = 0.,0.,0
    cΔy²,cR²         = maxΔy^2,maxR^2
    state            = allocate(pstate,Vector{State}(undef,saveiter ? maxiter : length(time))) # state is not a return argument so that data is not lost in case of exception
    for (step,t)     ∈ enumerate(time)
        verb && @printf "    step %3d" step
        old          = step==1 ? initial : state[step-1]
        s            = State(old.Λ,old.X,old.U,old.A,t,0.,model,dis)
        y            = s[dofgr] # includes scaling
        for iiter    = 1:maxiter
            citer   += 1
            asmt+=@elapsed assemble!(asm,dis,model,s, 0.,(dbg...,solver=:StaticX,step=step,iiter=iiter))
            solt+=@elapsed Δy = try asm.K\-asm.R catch; muscadeerror(@sprintf("Incremental solution failed at step=%i, iiter=%i",step,iiter)) end
            Δy²,R²   = sum(Δy.^2),sum(asm.R.^2)
            y      .+= Δy
            s[dofgr] = y  # includes descaling
            saveiter && (state[iiter]=s)
            if Δy²≤cΔy² && R²≤cR² 
                verb && @printf " converged in %3d iterations. |Δy|=%7.1e |R|=%7.1e\n" iiter √(Δy²) √(R²)
                saveiter || (state[step]=s)
                break#out of the iiter loop
            end
            iiter==maxiter && muscadeerror(@sprintf(" no convergence after %3d iterations |Δy|:%g / %g, |R|:%g / %g",iiter,√(Δy²),maxΔy,√(R²)^2,maxR))
        end
    end
    verb && @printf "\n    nel=%d, ndof=%d, nstep=%d, niter=%d, niter/nstep=%5.2f\n" getnele(model) getndof(dofgr) length(time) citer citer/length(time)
    verb && @printf "    Build  time = %s, (per iteration: %s, per iteration and element: %s)\n" showtime(asmt)  showtime(asmt/citer)  showtime(asmt/citer/getnele(model))
    verb && @printf "    Solve  time = %s, (per iteration: %s, per iteration and dof:     %s)\n" showtime(solt)  showtime(solt/citer)  showtime(solt/citer/getndof(dofgr))
end
