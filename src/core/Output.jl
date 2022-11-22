# TODO move to  modeldescription

function getdoftyp(model::Model,class::Symbol,field::Symbol)
    idoftyp = findfirst(doftyp.class==class && doftyp.field==field for doftyp∈model.doftyp)
    isnothing(idoftyp) && muscadeerror(@sprintf("The model has no dof of class %s and field %s.",class,field))    
    return model.doftyp[idoftyp]    
end
getdofID(model::Model,class::Symbol,field::Symbol) = getdoftyp(model,class,field).dofID
function getdofID(model::Model,class::Symbol,field::Symbol,nodID::AbstractVector{NodID})
    dofID  = getdofID(model,class,field) 
    i      = [model.dof[d].nodID ∈ nodID for d ∈ dofID]
    return dofID[i]
end

## Nodal results
function getdof_anc1(model,class,field,nodID,iders) 
    class ∈ [:Λ,:X,:U,:A] || muscadeerror(sprintf("Unknown dof class %s",class))
    c     = class==:Λ      ? :X                            : class
    dofID = nodID==NodID[] ? getdofID(model,c,field) : getdofID(model,c,field,nodID)
    iders = class∈[:Λ,:A]  ? [0]                           : iders
    return iders,dofID
end
function getdof_anc2!(out,state,class,iders,dofID)
    for ider∈iders
        s = if class==:Λ; state.Λ 
        elseif class==:X; state.X[ider+1]    
        elseif class==:U; state.U[ider+1]    
        elseif class==:A; state.A    
        end
        for (idof,d) ∈ enumerate(dofID)
            out[idof,ider+1] = s[d.idof] 
        end
    end 
end
function getdof(state::State;class::Symbol=:X,field::Symbol,nodID::Vector{NodID}=NodID[],iders::ℤ1=[0])  
    iders,dofID = getdof_anc1(state.model,class,field,nodID,iders)
    out         = Array{𝕣,2}(undef,length(dofID),length(iders)) # out[inod,ider+1]
    getdof_anc2!(out,state,class,iders,dofID)
    return out,dofID
end
function getdof(state::Vector{State};class::Symbol=:X,field::Symbol,nodID::Vector{NodID}=NodID[],iders::ℤ1=[0])
    iders,dofID = getdof_anc1(state[begin].model,class,field,nodID,iders)
    out         = Array{𝕣,3}(undef,length(dofID),length(iders),length(state)) # out[inod,ider+1]
    for istate ∈ eachindex(state)
        getdof_anc2!(view(out,:,:,istate),state[istate],class,iders,dofID)
    end
    return out,dofID
end





