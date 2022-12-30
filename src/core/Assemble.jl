
######## The disassembler

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
struct EletypDisassembler{nX,nU,nA}
    index :: Vector{XUA{𝕫,nX,nU,nA}}
    scale :: ΛXUA{𝕣,nX,nU,nA}
end
# dis.dis[ieletyp].index.[iele].X|U|A[ieledof]
# dis.dis[ieletyp].scale.Λ|X|U|A[ieledof]
# dis.scaleΛ|X|U|A[imoddof]
struct Disassembler{nX,nU,nA}
    dis     :: Vector{EletypDisassembler} 
    scaleΛ  :: 𝕣1
    scaleX  :: 𝕣1
    scaleU  :: 𝕣1
    scaleA  :: 𝕣1
end
function Disassembler(model::Model)
    neletyp                   = length(model.eleobj)  
    dis                       = Vector{EletypDisassembler}(undef,neletyp)
    nX,nU,nA                  = getndof(model,(:X,:U,:A))
    scaleΛ                    = Vector{𝕣}(undef,nX) # scale for state
    scaleX                    = Vector{𝕣}(undef,nX)
    scaleU                    = Vector{𝕣}(undef,nU)
    scaleA                    = Vector{𝕣}(undef,nA)
    for ieletyp               = 1:neletyp
        nele                  = length(model.eleobj[ieletyp])  
        E                     = eltype(model.eleobj[ieletyp])
        nX,nU,nA              = getndof(E,(:X,:U,:A))
        sΛ,sX,sU,sA           = 𝕣1(undef,nX),𝕣1(undef,nX),𝕣1(undef,nU),𝕣1(undef,nA) # scale for element
        ixdof,iudof,iadof     = 0,0,0
        for dofID             ∈ model.ele[ieletyp][begin].dofID
            doftyp            = getdoftyp(model,dofID)
            class,scale       = doftyp.class,doftyp.scale
            if     class == :X
                ixdof        += 1
                sX[ixdof]     = scale
                sΛ[ixdof]     = scale * model.scaleΛ
            elseif class == :U
                iudof        += 1
                sU[iudof]     = scale
            elseif class == :A
                iadof        += 1
                sA[iadof] = scale
            end
        end
        scale                 = ΛXUA{𝕣,nX,nU,nA}(sΛ,sX,sU,sA) # scale for element type
        iX,iU,iA              = 𝕫1(undef,nX),𝕫1(undef,nU),𝕫1(undef,nA)  # tmp arrays fof index into state of eledofs
        index                 = Vector{XUA{𝕫,nX,nU,nA}}(undef,nele)     # such indexes, for all elements in type
        for iele              = 1:nele
            ixdof,iudof,iadof = 0,0,0
            for dofID         ∈ model.ele[ieletyp][iele].dofID
                doftyp        = getdoftyp(model,dofID)
                class         = doftyp.class
                idof          = dofID.idof
                if     class == :X
                    ixdof    += 1
                    iX[ixdof] = idof  
                elseif class == :U
                    iudof    += 1
                    iU[iudof] = idof
                elseif class == :A
                    iadof    += 1
                    iA[iadof] = idof
                else
                    muscadeerror("element dof class must be :X,:U or :A")
                end
            end
            index[iele]       = XUA{𝕫,nX,nU,nA}(iX,iU,iA)
            scaleΛ[iX]        = scale.Λ  # "assemble" element scales into state scales
            scaleX[iX]        = scale.X
            scaleU[iU]        = scale.U
            scaleA[iA]        = scale.A
        end # for iele
        dis[ieletyp]          = EletypDisassembler{nX,nU,nA}(index,scale)
    end # for ieletyp
    return Disassembler(dis,scaleΛ,scaleX,scaleU,scaleA)
end

#### DofGroup

struct DofGroup{T1,T2,T3,T4,T5,T6,T7,T8} 
    nX     :: 𝕫 # of the _model_
    nU     :: 𝕫
    nA     :: 𝕣

    iΛ     :: T1   # state.Λ[iΛ] <-> y[jΛ]*Λscale
    iX     :: T2 
    iU     :: T3 
    iA     :: T4 

    jΛ     :: T5 
    jX     :: T6 
    jU     :: T7 
    jA     :: T8 

    Λscale :: 𝕣1
    Xscale :: 𝕣1
    Uscale :: 𝕣1
    Ascale :: 𝕣1
end
function DofGroup(dis::Disassembler,iΛ,iX,iU,iA) 
    # constructor for dofgroup with permutation within classe.  The datastructure of DofGroup supports dofgroups with arbitrary permutations - write another constructor
    nX,nU,nA    = length(dis.scaleX),length(dis.scaleU),length(dis.scaleA) # number of dofs in _model_
    nλ,nx,nu,na = length(iΛ),length(iX),length(iU),length(iA)              # number of dofs of each class in group
    jΛ,jX,jU,jA = gradientpartition(nλ,nx,nu,na)                               # we stack classes on top of each other in group vectors
    Λs,Xs,Us,As = dis.Λscale[iΛ],dis.Xscale[iX],dis.Uscale[iU],dis.Ascale[iA]
    return DofGroup(nX,nU,nA, iΛ,iX,iU,iA,  jΛ,jX,jU,jA, Λs,Xs,Us,As)
end
function decrement!(s::State,der::𝕫,y::𝕣1,gr::DofGroup) 
    for i ∈ eachindex(gr.iΛ); s.Λ[       gr.iΛ[i]] -= y[gr.jΛ[i]] * gr.Λscale[i]; end
    for i ∈ eachindex(gr.iX); s.X[der+1][gr.iX[i]] -= y[gr.jX[i]] * gr.Xscale[i]; end
    for i ∈ eachindex(gr.iU); s.U[der+1][gr.iU[i]] -= y[gr.jU[i]] * gr.Uscale[i]; end
    for i ∈ eachindex(gr.iA); s.A[       gr.iA[i]] -= y[gr.jA[i]] * gr.Ascale[i]; end
end
getndof(gr::DofGroup) = length(gr.iΛ)+length(gr.iX)+length(gr.iU)+length(gr.iA)
allΛdofs(model::Model,dis)   = DofGroup(dis, 1:getndof(model,:X),𝕫[],𝕫[],𝕫[])
allXdofs(model::Model,dis)   = DofGroup(dis, 𝕫[],1:getndof(model,:X),𝕫[],𝕫[])
allUdofs(model::Model,dis)   = DofGroup(dis, 𝕫[],𝕫[],1:getndof(model,:U),𝕫[])
allAdofs(model::Model,dis)   = DofGroup(dis, 𝕫[],𝕫[],𝕫[],1:getndof(model,:A))
allΛXUdofs(model::Model,dis) = DofGroup(dis, 1:getndof(model,:X),1:getndof(model,:X),1:getndof(model,:U),𝕫[])

# asm[iarray,ieletyp][idof/inz,iele] has value zero for gradient/hessian terms that are not to be added in.
#
# function prepare
#   allocate asm = Matrix{𝕫2}(undef,narray,neletyp)
#   for each array
#   pass @view asm[iarray,:] to preparevec/perparemat
# function preparevec/perparemat
#   for each ieletyp
#   asm[ieletyp] = Matrix{𝕫2}(undef,ndof/nnz,nele)
#   asm[ieletyp][idof/inz,iele] = ...
# function assemble!
#   for each ieletyp
#   pass @view asm[:,ieletyp] to assemblekernel!
# function assemblekernel!
#   for each iele
#   pass asm[:] and iele to addin!
# function addin!
#   for each array
#   use asm[iarray][:,iele] 

######## Prepare assemblers

function indexedstate(gr::DofGroup)
    # create a "state"  (Λ,X,U,A) of indices into the group - with zeros for modeldofs not in group
    Λ        = zeros(𝕫,gr.nX)
    X        = zeros(𝕫,gr.nX)
    U        = zeros(𝕫,gr.nU)
    A        = zeros(𝕫,gr.nA)
    Λ[gr.iΛ] = gr.jΛ
    X[gr.iX] = gr.jX
    U[gr.iU] = gr.jU
    A[gr.iA] = gr.jA
    return Λ,X,U,A
end
function gradientstructure(dofgr,dis::EletypDisassembler)
    # number of dofs of each class in the gradient returned by an element
    # because adiff is what it is, the gradient contains either all or no dofs in any given class
    nΛ       = length(dofgr.iΛ)==0 ? 0 : length(dis.scale.Λ) 
    nX       = length(dofgr.iX)==0 ? 0 : length(dis.scale.X) 
    nU       = length(dofgr.iU)==0 ? 0 : length(dis.scale.U) 
    nA       = length(dofgr.iA)==0 ? 0 : length(dis.scale.A) 
    return nΛ,nX,nU,nA
end
function gradientpartition(nΛ,nX,nU,nA)
    # indices into the class partitions of the gradient returned by an element
    iΛ          =          (1:nΛ)  
    iX          = nΛ+      (1:nX)
    iU          = nΛ+nX+   (1:nU)  
    iA          = nΛ+nX+nU+(1:nA)
    return iΛ,iX,iU,iA
end

# asm[ieletyp][idof|inz,iele] (its a @view)
# dofgr.iX[iXdof],dofgr.jX
# dis[ieletyp].index[iele].X|U|A[ieledof]
function preparevec!(asm,dofgr,dis) 
    # asm[ieletyp] == undef, please fill 
    Λ,X,U,A  = indexedstate(gr)                   # create a state of indices into the group - with zeros for modeldofs not in group
    for (ieletyp,di) ∈ enumerate(dis.dis)
        nΛ,nX,nU,nA = gradientstructure(dofgr,di) # number of dofs of each class in the gradient returned by an element
        iΛ,iX,iU,iU = gradientpartition(nΛ,nX,nU,nA)  # indices into said gradient
        # asm[ieletyp][idof,iele] (its a @view)
        asm[ieletyp] = zeros(𝕫,undef,nΛ+nX+nU+nA,length(di.index))
        for (iele,index) ∈ enumerate(di.index)
            asm[ieletyp][iΛ,iele] = Λ[index.X]
            asm[ieletyp][iX,iele] = X[index.X]
            asm[ieletyp][iU,iele] = U[index.U]
            asm[ieletyp][iA,iele] = A[index.A]
        end
    end
    return 𝕣1(undef,getndof(dofgr))
end
function preparemat!(asm,iasm,jasm,nidof,njdof) 
    # 1) traverse all eletyp
    #    compute number npair of contribution
    npair = 0
    for ieletyp ∈ eachindex(iasm)
        for iele = 1:size(iasm[ieletyp],2)
            npair += sum(iasm[ieletyp][:,iele].≠0)*sum(jasm[ieletyp][:,iele].≠0)
        end
    end
    # 2) traverse all elements 
    #       prepare a Vector A of all (jmoddof,imoddof) (in that order, for sort to work!) pairs of model dofs ::Vector{Tuple{Int64, Int64}}(undef,N)
    A = Vector{Tuple{𝕫,𝕫}}(undef,npair)
    ipair = 0
    for ieletyp ∈ eachindex(iasm)
        neledof,nele = size(iasm[ieletyp])
        for iele=1:nele, jeledof=1:neledof, ieledof=1:neledof
            if (iasm[ieletyp][ieledof,iele]≠0)  &&  (jasm[ieletyp][jeledof,iele]≠0)
                ipair += 1
                A[ipair] = (jasm[ieletyp][jeledof,iele] , iasm[ieletyp][ieledof,iele]) # NB: (j,i), not (i,j), because of lexicographic sortperm
            end
        end
    end
    # 3) sortperm(A)
    I = sortperm(A)
    # 4) traverse A[I] 
    #      find nnz
    #      create a list J that to each element of A associates an entry 1≤inz≤nnz into nzval
    #      prepare sparse
    J      = 𝕫1(undex,npair)
    nzval  = ones(𝕣,nnz) # could this be left undef?
    colptr = 𝕫1(undef,njdof+1)
    rowval = 𝕫1(undef,nnz)
    inz    = 0
    icol   = 1
    colptr[icol] = inz+1
    for ipair = 1:npair
        if (ipair==1) || (A[I[ipair]]≠A[I[ipair-1]]) 
            inz +=1
            (j,i) = A[I[ipair]] # NB: (j,i), not (i,j)
            rowval[inz] = i
            while j>icol
                icol +=1
                colptr[icol] = inz  # Column icol is in colptr[icol]:(colptr[icol+1]-1)
            end
        end
        J[ipair] = inz 
    end    
    # 5) traverse all elements again to distribute J into asm
    ipair = 0
    for ieletyp ∈ eachindex(iasm)
        neledof,nele = size(iasm[ieletyp])
        asm[ieletyp] = zeros(𝕫,neledof^2,nele)
        for iele=1:nele, jeledof=1:neledof, ieledof=1:neledof
            if (iasm[ieletyp][ieledof,iele]≠0)  &&  (jasm[ieletyp][jeledof,iele]≠0)
                ipair += 1
                ientry = ieledof+neledof*(jeledof-1) # TODO check transposition
                asm[ieletyp][ientry,iele] = J[ipair]  
            end
        end
    end
    # 6)
    return SparseMatrixCSC(nidof,njdof,colptr,rowval,nzval)   
end

####

function addinvec!(vec::Vector,asm,iele,array)
    for (i,eli) ∈ enumerate(array)
        j = asm[i,iele]
        if j≠0
            vec[j]+=eli
        end
    end
end   

######## Generic assembler

function assemble!(out,asm,dis,model,state,ε,dbg)
    zero!(out)
    for ieletyp ∈ eachindex(model.eleobj)
        eleobj  = model.eleobj[ieletyp]
        assemblesequential!(out,@view(asm,:,ieletyp),dis.dis[ieletyp], eleobj,state,ε,(dbg...,ieletyp=ieletyp))
    end
end
function assemblesequential!(out,asm,dis,eleobj,state,ε,dbg) 
    scale     = dis.scale
    for iele  ∈ eachindex(eleobj)
        index = dis.index[iele]
        Λe    = state.Λ[index.X]                 
        Xe    = Tuple(x[index.X] for x∈state.X)
        Ue    = Tuple(u[index.U] for u∈state.U)
        Ae    = state.A[index.A]
        addin!(out,asm,iele,scale,eleobj[iele],Λe,Xe,Ue,Ae, state.time,ε,(dbg...,iele=iele))
    end
end

###### scaled functions

function scaledlagrangian(scale,eleobj::E,Λs,Xs,Us,As, t,ε,dbg) where{E<:AbstractElement}
    Λ     =       Λs.*scale.Λ                 
    X     = Tuple(xs.*scale.X for xs∈Xs)  # TODO Tuple is slow, not typestable
    U     = Tuple(us.*scale.U for us∈Us)
    A     =       As.*scale.A
    L     = lagrangian(eleobj,Λ,X,U,A, t,ε,dbg)
    hasnan(L) && muscadeerror(dbg,"NaN in a Lagrangian or its partial derivatives")
    return L
end    
function scaledresidual(scale,eleobj::E, Xs,Us,As, t,ε,dbg) where{E<:AbstractElement} 
    X     = Tuple(xs.*scale.X for xs∈Xs)
    U     = Tuple(us.*scale.U for us∈Us)
    A     =       As.*scale.A
    R     = scale.Λ .* residual(eleobj, X,U,A, t,ε,dbg) 
    hasnan(R) && muscadeerror(dbg,"NaN in a residual or its partial derivatives")
    return R
end

