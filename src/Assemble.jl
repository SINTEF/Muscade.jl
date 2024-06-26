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
# dis.dis[ieletyp].index[iele].X|U|A[ieledof]       - disassembling model state into element dofs
# dis.dis[ieletyp].scale.Λ|X|U|A[ieledof]           - scaling each element type 
# dis.scaleΛ|X|U|A[imoddof]                         - scaling the model state
# dis.field  X|U|A[imoddof]                         - field of dofs in model state
struct Disassembler
    dis::Vector{EletypDisassembler}
    scaleΛ  :: 𝕣1
    scaleX  :: 𝕣1
    scaleU  :: 𝕣1
    scaleA  :: 𝕣1
    fieldX  :: Vector{Symbol}
    fieldU  :: Vector{Symbol}
    fieldA  :: Vector{Symbol}
end
function Disassembler(model::Model)
    neletyp                   = length(model.eleobj)  
    dis                       = Vector{EletypDisassembler}(undef,neletyp)
    NX,NU,NA                  = getndof(model,(:X,:U,:A))
    scaleΛ                    = Vector{𝕣}(undef,NX) # scale for state
    scaleX                    = Vector{𝕣}(undef,NX)
    scaleU                    = Vector{𝕣}(undef,NU)
    scaleA                    = Vector{𝕣}(undef,NA)
    fieldX                    = Vector{Symbol}(undef,NX)
    fieldU                    = Vector{Symbol}(undef,NU)
    fieldA                    = Vector{Symbol}(undef,NA)
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
                sA[iadof]     = scale
            end
        end
        scale                 = ΛXUA{𝕣,nX,nU,nA}(sΛ,sX,sU,sA) # scale for element type
        iX,iU,iA              = 𝕫1(undef,nX),𝕫1(undef,nU),𝕫1(undef,nA)  # tmp arrays for index into state of eledofs
        index                 = Vector{XUA{𝕫,nX,nU,nA}}(undef,nele)     # such indexes, for all elements in type
        for iele              = 1:nele
            ixdof,iudof,iadof = 0,0,0
            for dofID         ∈ model.ele[ieletyp][iele].dofID
                doftyp        = getdoftyp(model,dofID)
                class         = doftyp.class
                field         = doftyp.field
                idof          = dofID.idof  # model idof
                if     class == :X
                    ixdof    += 1
                    iX[ixdof] = idof  
                    fieldX[idof]= field
                elseif class == :U
                    iudof    += 1
                    iU[iudof] = idof
                    fieldU[idof]= field
                elseif class == :A
                    iadof    += 1
                    iA[iadof] = idof
                    fieldA[idof]= field
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
    return Disassembler(dis,scaleΛ,scaleX,scaleU,scaleA,fieldX,fieldU,fieldA)
end


######## state and initstate
# at each step, contains the complete, unscaled state of the system
mutable struct State{nΛder,nXder,nUder,TSP}
    time  :: 𝕣
    Λ     :: NTuple{nΛder,𝕣1}
    X     :: NTuple{nXder,𝕣1}
    U     :: NTuple{nUder,𝕣1}
    A     :: 𝕣1
    χ     :: Vector
    SP    :: TSP # solver parameter
    model :: Model
    dis   :: Disassembler
end
# a constructor that provides an initial state
State(model::Model,dis,χ;time=-∞) = State(time,(zeros(getndof(model,:X)),),
                                               (zeros(getndof(model,:X)),),
                                               (zeros(getndof(model,:U)),),
                                                zeros(getndof(model,:A))  ,
                                               χ,nothing,model,dis)
function State{nΛder,nXder,nUder}(s::State,SP::TSP) where{nΛder,nXder,nUder,TSP}
    Λ = ntuple(i->copy(∂n(s.Λ,i-1)),nΛder)
    X = ntuple(i->copy(∂n(s.X,i-1)),nXder)
    U = ntuple(i->copy(∂n(s.U,i-1)),nUder)
    State{nΛder,nXder,nUder,TSP}(s.time,Λ,X,U,copy(s.A),deepcopy(s.χ),SP,s.model,s.dis)
end 
State{nΛder,nXder,nUder}(s::State) where{nΛder,nXder,nUder} = State{nΛder,nXder,nUder}(s,(;))




#### DofGroup

# describes the relation between the dofs of the model, and a dof-vector containing an ordered selection
# of the dofs of the model.

struct DofGroup{T1,T2,T3,T4,T5,T6,T7,T8} 
    nX     :: 𝕫 # of the _model
    nU     :: 𝕫
    nA     :: 𝕫

    iΛ     :: T1   # state.Λ[iΛ] <-> y[jΛ]*scaleΛ (hence dofgroups can handle permutations)
    iX     :: T2 
    iU     :: T3 
    iA     :: T4 

    jΛ     :: T5 
    jX     :: T6 
    jU     :: T7 
    jA     :: T8 

    scaleΛ :: 𝕣1
    scaleX :: 𝕣1
    scaleU :: 𝕣1
    scaleA :: 𝕣1

    fieldΛ :: Vector{Symbol}  # fieldΛ[iΛ]
    fieldX :: Vector{Symbol}
    fieldU :: Vector{Symbol}
    fieldA :: Vector{Symbol}
end
function DofGroup(dis::Disassembler,iΛ,iX,iU,iA) 
    # constructor for dofgroup with permutation within each dof-class.  
    # The datastructure of DofGroup supports dofgroups with arbitrary permutations - but not this constructor
    nX,nU,nA    = length(dis.scaleX),length(dis.scaleU),length(dis.scaleA) # number of dofs in _model_
    nλ,nx,nu,na = length(iΛ),length(iX),length(iU),length(iA)              # number of dofs of each class in group
    jΛ,jX,jU,jA = gradientpartition(nλ,nx,nu,na)                           # we stack classes on top of each other in group vectors
    Λs,Xs,Us,As = dis.scaleΛ[iΛ],dis.scaleX[iX],dis.scaleU[iU],dis.scaleA[iA]
    Λf,Xf,Uf,Af = dis.fieldX[iΛ],dis.fieldX[iX],dis.fieldU[iU],dis.fieldA[iA]
    return DofGroup(nX,nU,nA, iΛ,iX,iU,iA,  jΛ,jX,jU,jA, Λs,Xs,Us,As, Λf,Xf,Uf,Af)
end
# use a dof-vector to decrement/increment/set/get the corresponding dofs in a State
function decrement!(s::State,der::𝕫,y::AbstractVector{𝕣},gr::DofGroup) 
    for i ∈ eachindex(gr.iΛ); s.Λ[der+1][gr.iΛ[i]] -= y[gr.jΛ[i]] * gr.scaleΛ[i]; end
    for i ∈ eachindex(gr.iX); s.X[der+1][gr.iX[i]] -= y[gr.jX[i]] * gr.scaleX[i]; end
    for i ∈ eachindex(gr.iU); s.U[der+1][gr.iU[i]] -= y[gr.jU[i]] * gr.scaleU[i]; end
    for i ∈ eachindex(gr.iA); s.A[       gr.iA[i]] -= y[gr.jA[i]] * gr.scaleA[i]; end
end
function increment!(s::State,der::𝕫,y::AbstractVector{𝕣},gr::DofGroup) 
    for i ∈ eachindex(gr.iΛ); s.Λ[der+1][gr.iΛ[i]] += y[gr.jΛ[i]] * gr.scaleΛ[i]; end
    for i ∈ eachindex(gr.iX); s.X[der+1][gr.iX[i]] += y[gr.jX[i]] * gr.scaleX[i]; end
    for i ∈ eachindex(gr.iU); s.U[der+1][gr.iU[i]] += y[gr.jU[i]] * gr.scaleU[i]; end
    for i ∈ eachindex(gr.iA); s.A[       gr.iA[i]] += y[gr.jA[i]] * gr.scaleA[i]; end
end
function set!(s::State,der::𝕫,y::AbstractVector{𝕣},gr::DofGroup) 
    s.Λ[der+1] .= 0
    s.X[der+1] .= 0
    s.U[der+1] .= 0
    s.A        .= 0
    for i ∈ eachindex(gr.iΛ); s.Λ[der+1][gr.iΛ[i]] = y[gr.jΛ[i]] * gr.scaleΛ[i]; end
    for i ∈ eachindex(gr.iX); s.X[der+1][gr.iX[i]] = y[gr.jX[i]] * gr.scaleX[i]; end
    for i ∈ eachindex(gr.iU); s.U[der+1][gr.iU[i]] = y[gr.jU[i]] * gr.scaleU[i]; end
    for i ∈ eachindex(gr.iA); s.A[       gr.iA[i]] = y[gr.jA[i]] * gr.scaleA[i]; end
end
function getdof!(s::State,der::𝕫,y::AbstractVector{𝕣},gr::DofGroup) 
    for i ∈ eachindex(gr.iΛ); y[gr.jΛ[i]] = s.Λ[der+1][gr.iΛ[i]] / gr.scaleΛ[i]; end
    for i ∈ eachindex(gr.iX); y[gr.jX[i]] = s.X[der+1][gr.iX[i]] / gr.scaleX[i]; end
    for i ∈ eachindex(gr.iU); y[gr.jU[i]] = s.U[der+1][gr.iU[i]] / gr.scaleU[i]; end
    for i ∈ eachindex(gr.iA); y[gr.jA[i]] = s.A[       gr.iA[i]] / gr.scaleA[i]; end
end
# create a tuple (Λ,X,U,A) of indices into the dofgroup - with zeros for modeldofs not in dofgroup
# so the model's iλ-th Λdof is found in y[Λ[iλ]]
function indexedstate(gr::DofGroup)
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
getndof(gr::DofGroup) = length(gr.iΛ)+length(gr.iX)+length(gr.iU)+length(gr.iA)

# some usefull Dofgroups
allΛdofs(  model::Model,dis) = DofGroup(dis, 1:getndof(model,:X),𝕫[],𝕫[],𝕫[])
allXdofs(  model::Model,dis) = DofGroup(dis, 𝕫[],1:getndof(model,:X),𝕫[],𝕫[])
allUdofs(  model::Model,dis) = DofGroup(dis, 𝕫[],𝕫[],1:getndof(model,:U),𝕫[])
allAdofs(  model::Model,dis) = DofGroup(dis, 𝕫[],𝕫[],𝕫[],1:getndof(model,:A))
allΛXUdofs(model::Model,dis) = DofGroup(dis, 1:getndof(model,:X),1:getndof(model,:X),1:getndof(model,:U),𝕫[])
allΛXUAdofs(model::Model,dis) = DofGroup(dis, 1:getndof(model,:X),1:getndof(model,:X),1:getndof(model,:U),1:getndof(model,:A))
function selecteddofs(model::Model,dis,classes)
    iΛ = :Λ ∈ classes ? (1:getndof(model,:X)) : 𝕫[] 
    iX = :X ∈ classes ? (1:getndof(model,:X)) : 𝕫[] 
    iU = :U ∈ classes ? (1:getndof(model,:U)) : 𝕫[] 
    iA = :A ∈ classes ? (1:getndof(model,:A)) : 𝕫[] 
    return DofGroup(dis, iΛ,iX,iU,iA)
end

######## Prepare assembler datastructure "asm"

# asm[iarray,ieletyp][ieledof/ientry,iele] has value zero for terms from element gradient/hessian that are not to be added in. Otherwise, the value they
# have is where in the matrix/vector/nzval to put the values.
# Example: for stiffness matrix iarray=2, beam element ieletyp=3, put the 4th entry (column major) of the iele=5th element into
# the asm[2,3][4,5]-th non-zero value (nzval) of the stiffness matrix for the solver.

# number of dofs of each class in the gradient returned by an element
# because adiff is what it is, the gradient contains either all or no dofs in any given class
function gradientstructure(dofgr,dis::EletypDisassembler)
    nΛ       = length(dofgr.iΛ)==0 ? 0 : length(dis.scale.Λ) 
    nX       = length(dofgr.iX)==0 ? 0 : length(dis.scale.X) 
    nU       = length(dofgr.iU)==0 ? 0 : length(dis.scale.U) 
    nA       = length(dofgr.iA)==0 ? 0 : length(dis.scale.A) 
    return nΛ,nX,nU,nA
end
# indices into the class partitions of the gradient returned by an element
function gradientpartition(nΛ,nX,nU,nA)
    iΛ          =           (1:nΛ)
    iX          = nΛ      .+(1:nX)
    iU          = nΛ+nX   .+(1:nU) 
    iA          = nΛ+nX+nU.+(1:nA)
    return iΛ,iX,iU,iA
end
# used in asmvec_kernel!
nonzeros(v) = v[v.≠0]
# prepare a vector-assembler.  Modifies asm.  Tyically called with view(asm,iarray,:), so inside the function
# asm is treated as asm[ieletyp][idof,iele]. Allocates and returns memory for the array to be assembled
function asmvec!(asm,dofgr,dis) 
    Λ,X,U,A  = indexedstate(dofgr)      # create a state of indices into the group - with zeros for modeldofs not in group
    for ieletyp ∈ eachindex(dis.dis)
        asmvec_kernel!(asm,ieletyp,dofgr,dis.dis[ieletyp],Λ,X,U,A)
    end
    return 𝕣1(undef,getndof(dofgr))
end
function asmvec_kernel!(asm,ieletyp,dofgr,dis,Λ,X,U,A) 
    nΛ,nX,nU,nA = gradientstructure(dofgr,dis) # number of dofs of each class in the gradient returned by an element
    iΛ,iX,iU,iA = gradientpartition(nΛ,nX,nU,nA)  # indices into said gradient TODO type unstable, barrier function
    asm[ieletyp] = zeros(𝕫,nΛ+nX+nU+nA,length(dis.index)) # asm[ieletyp][idof,iele] (its a view)
    for (iele,index) ∈ enumerate(dis.index)
        asm[ieletyp][iΛ,iele] .= nonzeros(Λ[index.X])  
        asm[ieletyp][iX,iele] .= nonzeros(X[index.X])
        asm[ieletyp][iU,iele] .= nonzeros(U[index.U])
        asm[ieletyp][iA,iele] .= nonzeros(A[index.A])
    end
end
function asmfullmat!(asm,iasm,jasm,nimoddof,njmoddof) 
    for ieletyp ∈ eachindex(iasm)
        nieledof,nele = size(iasm[ieletyp])
        njeledof      = size(jasm[ieletyp],1)
        asm[ieletyp]  = zeros(𝕫,nieledof*njeledof,nele)
        for iele=1:nele, jeledof=1:njeledof, ieledof=1:nieledof
            imoddof,jmoddof = iasm[ieletyp][ieledof,iele], jasm[ieletyp][jeledof,iele]
            if (imoddof≠0)  &&  (jmoddof≠0)
                ientry = ieledof+nieledof*(jeledof-1)
                asm[ieletyp][ientry,iele] = imoddof+nimoddof*(jmoddof-1)
            end
        end
    end
    return 𝕣2(undef,nimoddof,njmoddof)
end
function asmmat!(asm,iasm,jasm,nimoddof,njmoddof) 
    # 1) traverse all eletyp
    #    compute number npairs of contribution
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
        nieledof,nele = size(iasm[ieletyp])
        njeledof      = size(jasm[ieletyp],1)
        for iele=1:nele, jeledof=1:njeledof, ieledof=1:nieledof
            if (iasm[ieletyp][ieledof,iele]≠0)  &&  (jasm[ieletyp][jeledof,iele]≠0)
                ipair += 1
                A[ipair] = (jasm[ieletyp][jeledof,iele] , iasm[ieletyp][ieledof,iele]) # NB: (j,i), not (i,j), because of lexicographic sortperm
            end
        end
    end
    # 3) sortperm(A)
    I = sortperm(A)
    # 4) traverse A[I] 
    #      count nnz
    #      create a list J that to each element of A[I] associates an entry 1≤inz≤nnz into nzval
    #      create a list K that to each element of A    associates an entry 1≤inz≤nnz into nzval
    #      prepare sparse
    nnz    = 0
    for ipair = 1:npair
        if (ipair==1) || (A[I[ipair]]≠A[I[ipair-1]]) 
            nnz +=1
        end
    end    
    K      = 𝕫1(undef,npair) # to each pair in A    associate a unique entry number
    nzval  = ones(𝕣,nnz)     # could this be left undef and still get past the sparse constructor?  TODO 𝕣?
    colptr = 𝕫1(undef,njmoddof+1) # Column icol is in colptr[icol]:(colptr[icol+1]-1)
    colptr[1] = 1
    rowval = 𝕫1(undef,nnz)
    inz    = 0
    icol   = 1
    for ipair = 1:npair
        if (ipair==1) || (A[I[ipair]]≠A[I[ipair-1]]) 
            inz +=1
            (j,i) = A[I[ipair]] # NB: (j,i), not (i,j)
            rowval[inz] = i
            while j>icol 
                icol +=1
                colptr[icol] = inz  
            end
        end
        K[I[ipair]] = inz
    end    
    for i = icol+1:njmoddof+1
        colptr[i] = nnz+1
    end
    
    # 5) traverse all elements again to distribute J into asm
    ipair = 0
    for ieletyp ∈ eachindex(iasm)
        nieledof,nele = size(iasm[ieletyp]  )
        njeledof      = size(jasm[ieletyp],1)
        asm[ieletyp]  = zeros(𝕫,nieledof*njeledof,nele)
        for iele=1:nele, jeledof=1:njeledof, ieledof=1:nieledof
            if (iasm[ieletyp][ieledof,iele]≠0)  &&  (jasm[ieletyp][jeledof,iele]≠0)
                ipair += 1
                ientry = ieledof+nieledof*(jeledof-1) 
                asm[ieletyp][ientry,iele] = K[ipair]  
            end
        end
    end
    # 6)
    return SparseMatrixCSC(nimoddof,njmoddof,colptr,rowval,nzval)   
end


######## Assembler methods
# The call stack     who dunnit               why
# 
# MySolver           Muscade/SomeSolver.jl    wants some matrices and vectors
# assemble!          Muscade/Assemble.jl      loop over element types (barrier function)
# assemble_!         Muscade/Assemble.jl      loop over elements within type (typestable)
# addin!             Muscade/SomeSolver.jl    do adiff and add-in, clean χ from adiff, solver-specific
# getresidual        Muscade/Assemble.jl      typechecking the call from addin! call residual or Lagrangian as available, check for NaNs
# residual           MyElement.jl             the element code 


abstract type Assembly end # solver define concrete "assemblies" which is a collection of matrices and solvers wanted for a phase in the solution process


# sequential, called by the solver
function assemble!(out::Assembly,asm,dis,model,state,dbg) 
    zero!(out)
    for ieletyp = 1:lastindex(model.eleobj)
        eleobj  = model.eleobj[ieletyp]
        assemble_!(out,view(asm,:,ieletyp),dis.dis[ieletyp],eleobj,state,state.SP,(dbg...,ieletyp=ieletyp))
    end
end
function assemble_!(out::Assembly,asm,dis,eleobj,state::State{nΛder,nXder,nUder},SP,dbg) where{nΛder,nXder,nUder}
    scale     = dis.scale
    for iele  = 1:lastindex(eleobj)
        index = dis.index[iele]
        Λe    = NTuple{nΛder}(λ[index.X] for λ∈state.Λ)
        Xe    = NTuple{nXder}(x[index.X] for x∈state.X)
        Ue    = NTuple{nUder}(u[index.U] for u∈state.U)
        Ae    = state.A[index.A]
        addin!(out,asm,iele,scale,eleobj[iele],Λe,Xe,Ue,Ae, state.time,SP,(dbg...,iele=iele)) # defined by solver.  Called for each element. But the asm that is passed
    end                                                                                       # is of the form asm[iarray][ientry,iele], because addin! will add to all arrays in one pass
end

# multithreaded
using Base.Threads
one_for_each_thread(x) = SVector{nthreads()}(deepcopy(x) for i=1:nthreads())
firstelement(x::AbstractVector) = x[1]
firstelement(x                ) = x

function add!(a::Array,b::Array)
    for i∈eachindex(a)
        a[i] += b[i]
    end
end
function add!(a::SparseMatrixCSC,b::SparseMatrixCSC) # assumes identical sparsity structure
    for i∈eachindex(a.nzval)
        a.nzval[i] += b.nzval[i]
    end
end

function assemble!(out::AbstractVector{A},asm,dis,model,state,dbg) where{A<:Assembly}
    for i = 1:nthreads() 
        zero!(out[i])
    end
    for ieletyp = 1:lastindex(model.eleobj)
        eleobj  = model.eleobj[ieletyp]
        assemble_!(out,view(asm,:,ieletyp),dis.dis[ieletyp], eleobj,state,state.SP,(dbg...,ieletyp=ieletyp))
    end
    for i = 2:nthreads() 
        add!(out[1],out[i])
    end
end

function assemble_!(out::AbstractVector{A},asm,dis,eleobj,state::State{nΛder,nXder,nUder},SP,dbg) where{nΛder,nXder,nUder,A<:Assembly}
    scale     = dis.scale
    @threads for iele  = 1:lastindex(eleobj)
        index = dis.index[iele]
        Λe    = NTuple{nΛder}(λ[index.X] for λ∈state.Λ)
        Xe    = NTuple{nXder}(x[index.X] for x∈state.X)
        Ue    = NTuple{nUder}(u[index.U] for u∈state.U)
        Ae    = state.A[index.A]
        addin!(out[threadid()],asm,iele,scale,eleobj[iele],Λe,Xe,Ue,Ae, state.time,SP,(dbg...,iele=iele))
    end
end


############# Tools for addin!



#### zero!
"""
    zero!(a)

Set to zero all elements of an arrays. If `a` is sparse, 
the vector `nzval` of values is set to zero and the sparsity structure is unchanged.
""" 

function zero!(out::DenseArray)
    for i∈eachindex(out)
        out[i] = 0
    end
end
function zero!(out::AbstractSparseArray)
    for i∈eachindex(out.nzval)
        out.nzval[i] = 0
    end
end

#### extract value or derivatives from a SVector 'a' of adiffs, and add it directly into vector, full matrix or nzval of sparse matrix 'out'.
function add_value!(out::𝕣1,asm,iele,a::SVector{M,∂ℝ{P,N,𝕣}},ias) where{P,N,M}
    # asm[ientry,iel]
    for (ientry,ia) ∈ enumerate(ias)
        iout = asm[ientry,iele]
        if iout≠0
            out[iout]+=a[ia].x
        end
    end
end   
function add_value!(out::𝕣1,asm,iele,a::SVector{M,𝕣},ias) where{M}
    for (ientry,ia) ∈ enumerate(ias)
        iout = asm[ientry,iele]
        if iout≠0
            out[iout]+=a[ia]
        end
    end
end   
add_value!(out,asm,iele,a) = add_value!(out,asm,iele,a,eachindex(a)) 
struct add_∂!{P} end 
function add_∂!{P}(out::Array,asm,iele,a::SVector{M,∂ℝ{P,N,R}},i1as,i2as) where{P,N,R,M}
    for (i1asm,i1a) ∈ enumerate(i1as), (i2asm,i2a) ∈ enumerate(i2as)
        ientry = i1asm+length(i1as)*(i2asm-1)
        iout = asm[ientry,iele]
        if iout≠0
            out[iout]+=a[i1a].dx[i2a]  
        end
    end
end  
add_∂!{P}(out::SparseMatrixCSC,args...) where{P}                      = add_∂!{P}(out.nzval,args...)
add_∂!{P}(out::Array,asm,iele,a::SVector{M,R},args...) where{P,M,R}   = nothing
add_∂!{P}(out::Array,asm,iele,a::SVector{M,∂ℝ{P,N,R}}) where{P,N,R,M} = add_∂!{P}(out,asm,iele,a,SVector{M}(1:M),SVector{N}(1:N))


####### called by addin!, and by nested elements to "get a Lagrangian" and "get a residual"
# 1) comprehensive check of the types of arguments, to help catch bugs in solvers and elements
# 2) if a residual is wanted by the solver and only Lagrangian is implemented by the element (or the other way around), do some magic
# 3) check for NaNs in the results 
#
# Note that getLagrangian receives Λ::SVector. addin! by contrast receives Λ::NTuple{SVector}, this is not a bug

function getresidual(eleobj::Eleobj,  
    X::NTuple{Ndx,SVector{Nx,Rx}},
    U::NTuple{Ndu,SVector{Nu,Ru}},
    A::           SVector{Na,Ra} ,
    t::ℝ,χ,SP,dbg,req...)     where{Eleobj<:AbstractElement,Ndx,Nx,Rx<:ℝ,Ndu,Nu,Ru<:ℝ,Na,Ra<:ℝ} 

    if hasmethod(residual  ,(Eleobj,       NTuple,NTuple,𝕣1,𝕣,Any,Function,NamedTuple,NamedTuple))
        R,χn,FB,eleres... = residual(  eleobj,  X,U,A,t,χ,SP,dbg,req...)
        hasnan(R,χn,FB) && muscadeerror((dbg...,t=t,SP=SP),@sprintf("residual(%s,...) returned NaN in R, χ, FB or derivatives",Eleobj))  

    elseif hasmethod(lagrangian,(Eleobj,NTuple,NTuple,NTuple,𝕣1,𝕣,Any,Function,NamedTuple,NamedTuple))
        P   = constants(∂0(X),∂0(U),A,t)
        Λ   = δ{P,Nx,𝕣}() 
        L,χn,FB,eleres... = lagrangian(eleobj,Λ,X,U,A,t,χ,SP,dbg,req...)    
        hasnan(L,χn,FB) && muscadeerror((dbg...,t=t,SP=SP),@sprintf("lagrangian(%s,...) returned NaN in L, χ, FB or derivatives",Eleobj))   
        R = ∂{P,Nx}(L)
    else 
        muscadeerror((dbg...,t=t,SP=SP),@sprintf("Element %s must have method 'Muscade.lagrangian' or/and 'Muscade.residual' with correct interface",Eleobj))
    end
    return R,χn,FB,eleres...
end


function getlagrangian(eleobj::Eleobj,  
    Λ::           SVector{Nx,Rλ} ,  
    X::NTuple{Ndx,SVector{Nx,Rx}},
    U::NTuple{Ndu,SVector{Nu,Ru}},
    A::           SVector{Na,Ra} ,
    t::ℝ,χ,SP,dbg,req...)     where{Eleobj<:AbstractElement,Rλ<:ℝ,Ndx,Nx,Rx<:ℝ,Ndu,Nu,Ru<:ℝ,Na,Ra<:ℝ} 

    if     hasmethod(lagrangian,(Eleobj,NTuple,NTuple,NTuple,𝕣1,𝕣,Any,Function,NamedTuple,NamedTuple))
        L,χn,FB,eleres... = lagrangian(eleobj,Λ,X,U,A,t,χ,SP,dbg,req...)
        hasnan(L,χn,FB) && muscadeerror((dbg...,t=t,SP=SP),@sprintf("lagrangian(%s,...) returned NaN in L, χ, FB or derivatives",Eleobj))     
    elseif hasmethod(residual  ,(Eleobj,       NTuple,NTuple,𝕣1,𝕣,Any,Function,NamedTuple,NamedTuple))
        R,χn,FB,eleres... = residual(  eleobj,  X,U,A,t,χ,SP,dbg,req...)
        hasnan(R,χn,FB) && muscadeerror((dbg...,t=t,SP=SP),@sprintf("residual(%s,...) returned NaN in R, χ, FB or derivatives",Eleobj)) 
        L = Λ ∘₁ R
    else
        muscadeerror((dbg...,t=t,SP=SP),@sprintf("Element %s must have method 'Muscade.lagrangian' or/and 'Muscade.residual' with correct interface",Eleobj))
    end
    return L,χn,FB,eleres... 
end
