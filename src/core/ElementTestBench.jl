using  Printf,ForwardDiff,StaticArrays
using  Muscade # to define dofid


function testStaticElement(el; δX,X,U,A, t::Float64=0.,ε::Float64=0., verbose::Bool=true)
    id           = dofid(el)
#    n            = neldof(el) 
    dbg          = ()
    nX,nU,nA     = length.((X,U,A))
    iδX,iX,iU,iA = (1:nX) , (1:nX) .+ nX , (1:nU) .+ 2nX , (1:nA) .+ (2nX+nU)         
    closure(Y)   = lagrangian(el,Y[iδX],[Y[iX]],[Y[iU]],Y[iA], t,ε,dbg)
    Ly           = ForwardDiff.gradient(closure,vcat(δX,X,U,A))
    Lδx,Lx,Lu,La = Ly[iδX], Ly[iX], Ly[iU], Ly[iA]
    if verbose
        @printf "\nElement type: %s\n" typeof(el)
        if n.X > 0
            @printf "\n    idof               doftyp   inod          δX           X         Lδx          Lx \n"
            for idof = 1:n.X
                @printf "    %4d     %16s  %5d  %10.3g  %10.3g  %10.3g  %10.3g\n" idof id.X.typ[idof] id.X.nod[idof] δX[idof] X[idof] Lδx[idof] Lx[idof]
            end
        end
        if n.U > 0
            @printf "\n    idof               doftyp   inod           U          Lu \n"
            for idof = 1:n.U
                @printf "    %4d     %16s  %5d  %10.3g  %10.3g\n" idof id.U.typ[idof] id.U.nod[idof] U[idof] Lu[idof]
            end
        end
        if n.A > 0
            @printf "\n    idof               doftyp   inod           A          La \n"
            for idof = 1:n.A
                @printf "    %4d     %16s  %5d  %10.3g  %10.3g\n" idof id.A.typ[idof] id.A.nod[idof] A[idof] La[idof]
            end
        end
    end

    return Lδx,Lx,Lu,La
end


