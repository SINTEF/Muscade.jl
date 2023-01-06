## Nodal results
function getdof(state::State;kwargs...)  
    out,dofID = getdof([state];kwargs...)
    return reshape(out,size(out)[1:2]),dofID
end
function getdof(state::Vector{S};class::Symbol=:X,field::Symbol,nodID::Vector{NodID}=NodID[],iders::ℤ1=[0])where {S<:State}
    class ∈ [:Λ,:X,:U,:A] || muscadeerror(sprintf("Unknown dof class %s",class))
    c     = class==:Λ      ? :X                                   : class
    dofID = nodID==NodID[] ? getdofID(state[begin].model,c,field) : getdofID(state[begin].model,c,field,nodID)
    iders = class∈[:Λ,:A]  ? [0]                                  : iders
    out   = Array{𝕣,3}(undef,length(dofID),length(iders),length(state)) # out[inod,ider+1]
    for istate ∈ eachindex(state)
        for ider∈iders
            s = if class==:Λ; state[istate].Λ 
            elseif class==:X; state[istate].X[ider+1]    
            elseif class==:U; state[istate].U[ider+1]    
            elseif class==:A; state[istate].A    
            end
            for (idof,d) ∈ enumerate(dofID)
                out[idof,ider+1,istate] = s[d.idof] 
            end
        end 
    end
    return out,dofID
end

# Elemental results
function extractkernel!(out,key,eleobj,eleID,dis::EletypDisassembler,state::State,dbg) # typestable kernel
    for (iele,ei) ∈ enumerate(eleID)
        index = dis.index[ei.iele]
        Λ     = state.Λ[index.X]                 
        X     = Tuple(x[index.X] for x∈state.X)
        U     = Tuple(u[index.U] for u∈state.U)
        A     = state.A[index.A]
        _     = lagrangian(view(out,:,iele),key,eleobj[ei.iele],Λ,X,U,A,state.time,state.ε,(dbg...,iele=ei.iele))
    end
end
function getresult(state::Vector{S},req; eleID::Vector{EleID})where {S<:State}
    # One element type, some or all elements within the types
    # out[ikey,iele,istep]
    ieletyp             = eleID[begin].ieletyp
    all(e.ieletyp== ieletyp for e∈eleID) || muscadeerror("All elements must be of the same element type")
    eleobj              = state[begin].model.eleobj[ieletyp]
    dis                 = state[begin].dis.dis[ieletyp]
    key,nkey            = makekey(req,espyable(eltype(eleobj)))
    nstep,nele          = length(state),length(eleID)
    out                 = Array{𝕣,3}(undef,nkey,nele,nstep)
    for (istep,s) ∈ enumerate(state)
        extractkernel!(view(out,:,:,istep),key,eleobj,eleID,dis,s,(ieletyp=ieletyp,istep=istep))
    end
    return out,key
end
function getresult(state::State,req;kwargs...)  
    out,key = getresult([state],req;kwargs...)
    return reshape(out,size(out)[1:2]),key
end

function describeX(state::State)
    model = state.model
    nX    = getndof(model,:X)
    nder  = length(state.X)
    for iX = 1:nX
        dofID   = DofID(:X,iX)
        dof     = model.dof[dofID] 
        doftyp  = model.doftyp[dof.idoftyp]
        @printf "NodID(%i), class=:%s, field=:%-15s   " dof.nodID.inod dofID.class doftyp.field
        for ider = 1:nder
            @printf "%15g " state.X[ider][iX]
        end
        @printf "\n" 
    end
end
function describeΛX(state::State)
    model = state.model
    nX    = getndof(model,:X)
    nder  = length(state.X)
    for iX = 1:nX
        dofID   = DofID(:X,iX)
        dof     = model.dof[dofID] 
        doftyp  = model.doftyp[dof.idoftyp]
        @printf "NodID(%i), class=:%s, field=:%-15s   %15g " dof.nodID.inod dofID.class doftyp.field state.Λ[iX]
        for ider = 1:nder
            @printf "%15g " state.X[ider][iX]
        end
        @printf "\n" 
    end
end
function describeU(state::State)
    model = state.model
    nU    = getndof(model,:U)
    nder  = length(state.U)
    for iU = 1:nU
        dofID   = DofID(:U,iU)
        dof     = model.dof[dofID] 
        doftyp  = model.doftyp[dof.idoftyp]
        @printf "NodID(%i), class=:%s, field=:%-15s   " dof.nodID.inod dofID.class doftyp.field
        for ider = 1:nder
            @printf "%15g " state.U[ider][iU]
        end
        @printf "\n"
    end
end
function describeA(state::State)
    model = state.model
    nA    = getndof(model,:A)
    for iA = 1:nA
        dofID   = DofID(:A,iA)
        dof     = model.dof[dofID] 
        doftyp  = model.doftyp[dof.idoftyp]
        @printf "NodID(%i), class=:%s, field=:%-15s   %15g\n" dof.nodID.inod dofID.class doftyp.field state.A[iA] 
    end
end
function describe(state::State;class::Symbol=:all)
    if class ==:all
        describeΛX(state)
        describeU(state)
        describeA(state)
    elseif class==:Λ || class==:ΛX
        describeΛX(state)
    elseif class ==:X    
        describeX(state)
    elseif class ==:U    
        describeU(state)
    elseif class ==:A    
        describeA(state)
    else
        printstyled("Not a valid class\n",color=:red,bold=true)
    end
end