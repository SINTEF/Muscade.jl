module FloaterMotions

using StaticArrays, LinearAlgebra
using Muscade
export FloaterOnCalmWater

const rigid  = (:surge,:sway,:heave,:roll,:pitch,:yaw)
const idx    = (:11,:12,:13,:14,:15,:16,:22,:23,:24,:25,:26,:33,:34,:35,:36,:44,:45,:46,:55,:56,:66)
fold(x::SVector{21}) = SMatrix{6,6}(x[ 1],x[ 2],x[ 3],x[ 4],x[ 5],x[ 6],
                                    x[ 2],x[ 7],x[ 8],x[ 9],x[10],x[11],
                                    x[ 3],x[ 8],x[12],x[13],x[14],x[15],
                                    x[ 4],x[ 9],x[13],x[16],x[17],x[18],
                                    x[ 5],x[10],x[14],x[17],x[19],x[20],
                                    x[ 6],x[11],x[15],x[18],x[20],x[21])

struct FloaterOnCalmWater <: AbstractElement
    K  :: 𝕣2
    Cₒ :: 𝕣2
    Mₒ :: 𝕣2
end
FloaterOnCalmWater(nod::Vector{Node};K,Cₒ,Mₒ  )  = FloaterOnCalmWater(K,Cₒ,Mₒ  )

Muscade.doflist(::Type{<:FloaterOnCalmWater}) = (inod = (ntuple(i-> 1,6)...,ntuple(i-> 1,6)...,ntuple(i-> 1,21)...               ,ntuple(i-> 1,21)...               )                                             
                                                 class= (ntuple(i->:X,6)...,ntuple(i->:U,6)...,ntuple(i->:A,21)...               ,ntuple(i->:A,21)...               ), 
                                                 field= (rigid...          ,rigid...          ,ntuple(i->Symbol(:C,idx[i]),21)...,ntuple(i->Symbol(:M,idx[i]),21)...))

@espy function Muscade.residual(o::FloaterOnCalmWater,   X,U,A,t,SP,dbg) 
    x,x′,x″    = ∂0(X),∂1(X),∂2(X)   
    ☼u         = ∂0(U)
    a          = exp10.(A)
    ☼r₀        = o.K∘x
    ☼r₁        = (o.C₀.*fold(a[ 1:21]))∘x′
    ☼r₂        = (o.M₀.*fold(a[22:42]))∘x″
    return r₀+r₁+r₂-u,  noFB
end

end