using  Printf
using  Muscade 


function test_static_element(ele::eletyp; δX,X,U,A, t::Float64=0.,ε::Float64=0., verbose::Bool=true,dbg = ()) where{eletyp<:AbstractElement}
    inod,class,field = Muscade.getdoflist(eletyp)
    iXdof            = Muscade.getidof(eletyp,:X)
    iUdof            = Muscade.getidof(eletyp,:U)
    iAdof            = Muscade.getidof(eletyp,:A)
    nX,nU,nA         = Muscade.getndofs(eletyp)
    L,Lδx,Lx,Lu,La   = Muscade.gradient(SeverΛXUAstatic,ele,δX,[X],[U],A, t,ε,dbg)

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


