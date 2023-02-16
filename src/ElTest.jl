using  Printf
using  Muscade 

####### For testing: get all the gradients. 
function gradient(eleobj,Λ,X,U,A, t,γ,dbg) 
    P            = constants(Λ,∂0(X),∂0(U),A,t)
    nX,nU,nA     = length(Λ),length(∂0(U)),length(A)
    N            = 2nX+nU+nA
    iΛ,iX,iU,iA  = (1:nX) , (1:nX) .+ nX , (1:nU) .+ 2nX , (1:nA) .+ (2nX+nU)  
    ΔY           = δ{P,N,𝕣}()                        
    L,minγfac    = Muscade.getlagrangian(Muscade.implemented(eleobj)...,eleobj,Λ+ΔY[iΛ],(∂0(X)+ΔY[iX],),(∂0(U)+ΔY[iU],),A+ΔY[iA], t,γ,dbg)
    Ly           = ∂{P,N}(L)
    return (L=value{P}(L), Lλ=Ly[iΛ], Lx=Ly[iX], Lu=Ly[iU], La=Ly[iA])
end

function test_static_element(ele::eletyp; δX,X,U,A, t::Float64=0.,γ::Float64=0., verbose::Bool=true,dbg = NamedTuple()) where{eletyp<:AbstractElement}
    inod,class,field = Muscade.getdoflist(eletyp)
    iXdof            = Muscade.getidof(eletyp,:X)
    iUdof            = Muscade.getidof(eletyp,:U)
    iAdof            = Muscade.getidof(eletyp,:A)
    nX,nU,nA         = Muscade.getndof(eletyp,(:X,:U,:A))
    L,Lδx,Lx,Lu,La   = gradient(ele,δX,[X],[U],A, t,γ,dbg)

    if verbose
        @printf "\nElement type: %s\n" typeof(el)
        if nX > 0
            @printf "\n    idof               doftyp   inod          δX           X         Lδx          Lx \n"
            for iX = 1:nX
                idof = iXdof[iX]
                @printf "    %4d     %16s  %5d  %10.3g  %10.3g  %10.3g  %10.3g\n" idof field[idof] inod[idof] δX[idof] X[idof] Lδx[idof] Lx[idof]
            end
        end
        if nU > 0
            @printf "\n    idof               doftyp   inod           U          Lu \n"
            for iU = 1:nU
                idof = iUdof[iU]
                @printf "    %4d     %16s  %5d  %10.3g  %10.3g\n" idof field[idof] inod[idof] U[idof] Lu[idof]
            end
        end
        if nA > 0
            @printf "\n    idof               doftyp   inod           A          La \n"
            for iA = 1:nA
                idof = iAdof[iA]
                @printf "    %4d     %16s  %5d  %10.3g  %10.3g\n" idof field[idof] inod[idof] A[idof] La[idof]
            end
        end
    end

    return Lδx,Lx,Lu,La
end


