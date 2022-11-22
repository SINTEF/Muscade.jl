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
function getdof(state::State;kwargs...)  
    out,dofID = getdof([state];kwargs...)
    return reshape(out,size(out)[1:2]),dofID
end
function getdof(state::Vector{S};class::Symbol=:X,field::Symbol,nodID::Vector{NodID}=NodID[],iders::ℤ1=[0])where {S<:State}
    class ∈ [:Λ,:X,:U,:A] || muscadeerror(sprintf("Unknown dof class %s",class))
    c     = class==:Λ      ? :X                            : class
    dofID = nodID==NodID[] ? getdofID(state[begin].model,c,field) : getdofID(state[begin].model,c,field,nodID)
    iders = class∈[:Λ,:A]  ? [0]                           : iders
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





