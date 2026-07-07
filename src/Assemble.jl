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
    SP    :: TSP # solver parameter
    model :: Model
    dis   :: Disassembler
    # Inner constructors
    # Provide values, infer type
    State(time::𝕣, Λ::NTuple{nΛder,𝕣1}, X::NTuple{nXder,𝕣1}, U::NTuple{nUder,𝕣1}, A::𝕣1, SP::TSP, model::Model, dis::Disassembler) where{nΛder,nXder,nUder,TSP} =
        new{nΛder,nXder,nUder,TSP}(time,Λ,X,U,A,SP,model,dis) 
    # Provide type, undef'd values    
    State{nΛder,nXder,nUder,TSP}() where{nΛder,nXder,nUder,TSP} = new{nΛder,nXder,nUder,TSP}()   
end

# a constructor that provides an initial zero state, specify derivatives
State{nΛder,nXder,nUder}(model::Model,dis;time=-∞) where{nΛder,nXder,nUder} = 
                                  State(time,ntuple(i->zeros(getndof(model,:X)),nΛder),
                                             ntuple(i->zeros(getndof(model,:X)),nXder),
                                             ntuple(i->zeros(getndof(model,:U)),nUder),
                                                       zeros(getndof(model,:A))       ,
                                             nothing,model,dis)
# a shallow copy "constructor" to shave off unwanted derivatives (or pad with zeros) 
State{nΛder,nXder,nUder}(s::State,SP=s.SP) where{nΛder,nXder,nUder} = State{nΛder,nXder,nUder}(s.time,s.Λ,s.X,s.U,s.A,SP,s.model,s.dis)
# the same but from components
function State{nΛder,nXder,nUder}(time,Λ,X,U,A,SP::TSP,model,dis) where{nΛder,nXder,nUder,TSP}
    state       = State{nΛder,nXder,nUder,TSP}()
    state.time  = time
    state.Λ     = ntuple(i->∂n(Λ,i-1),nΛder)
    state.X     = ntuple(i->∂n(X,i-1),nXder)
    state.U     = ntuple(i->∂n(U,i-1),nUder)
    state.A     = A
    state.SP    = SP
    state.model = model
    state.dis   = dis
    return state
end
# the same but from components
function State{nΛder,nXder,nUder}(state;time=state.time,Λ=state.Λ,X=state.X,U=state.U,A=state.A,SP::TSP=state.SP,model=state.model,dis=state.dis) where{nΛder,nXder,nUder,TSP}
    state       = State{nΛder,nXder,nUder,TSP}()
    state.time  = time
    state.Λ     = ntuple(i->∂n(Λ,i-1),nΛder)
    state.X     = ntuple(i->∂n(X,i-1),nXder)
    state.U     = ntuple(i->∂n(U,i-1),nUder)
    state.A     = A
    state.SP    = SP
    state.model = model
    state.dis   = dis
    return state
end
"""
    s2 = copy(s1::State)

Copy a state so that `s1` and `s2` do not share memory.

(Actualy the `State`s will stil share the same `Model`, but no operations are provided to mutate it).

    s2 = copy(s1;[time=t])

create a copied `State` with time `t`. 

See also: [`initialize!`](@ref), [`setdof!`](@ref)
"""
Base.copy(s::State;time=s.time,SP=s.SP) = State(time,deepcopy(s.Λ),deepcopy(s.X),deepcopy(s.U),deepcopy(s.A),SP,s.model,s.dis) 

function Base.show(io::IO,o::State{nΛder,nXder,nUder,TSP}) where{nΛder,nXder,nUder,TSP}
    print(io,"State{",nΛder,",",nXder,",",nUder,",",TSP,"}(model=:",o.model.ID,", time=",o.time,", ...)")
end
function Base.show(io::IO,o::Vector{State{nΛder,nXder,nUder,TSP}}) where{nΛder,nXder,nUder,TSP}
    println(io,length(o),"-element Vector{State{",nΛder,",",nXder,",",nUder,",",TSP,"}}:")
    for (i,oᵢ) ∈ enumerate(o)
        println(io," ",i,": model=:",oᵢ.model.ID,", time=",oᵢ.time)
    end
end




#### DofGroup

# describes the relation between the dofs of the model, and a dof-vector containing an ordered selection
# of the dofs of the model.

struct DofGroup 
    nX     :: 𝕫 # of the _model_
    nU     :: 𝕫
    nA     :: 𝕫

    iΛ     :: 𝕫1   # state.Λ[iΛ] <-> y[jΛ]*scaleΛ (hence dofgroups can handle permutations)
    iX     :: 𝕫1 
    iU     :: 𝕫1 
    iA     :: 𝕫1 

    jΛ     :: 𝕫1 
    jX     :: 𝕫1 
    jU     :: 𝕫1 
    jA     :: 𝕫1 

    scaleΛ :: 𝕣1              # scaleΛ[iΛ]
    scaleX :: 𝕣1
    scaleU :: 𝕣1
    scaleA :: 𝕣1

    fieldΛ :: Vector{Symbol}  # fieldΛ[iΛ]
    fieldX :: Vector{Symbol}
    fieldU :: Vector{Symbol}
    fieldA :: Vector{Symbol}
    DofGroup(nX,nU,nA, iΛ,iX,iU,iA,  jΛ,jX,jU,jA, Λs,Xs,Us,As, Λf,Xf,Uf,Af) = 
      new(nX,nU,nA, collect(iΛ),collect(iX),collect(iU),collect(iA),  collect(jΛ),collect(jX),collect(jU),collect(jA), Λs,Xs,Us,As, Λf,Xf,Uf,Af)
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
function decrement!(s::State,ider::𝕫,y::AbstractVector{𝕣},gr::DofGroup) 
    if ider≤length(s.Λ) for i ∈ eachindex(gr.iΛ); s.Λ[ider][gr.iΛ[i]] -= y[gr.jΛ[i]] * gr.scaleΛ[i]; end end
    if ider≤length(s.X) for i ∈ eachindex(gr.iX); s.X[ider][gr.iX[i]] -= y[gr.jX[i]] * gr.scaleX[i]; end end
    if ider≤length(s.U) for i ∈ eachindex(gr.iU); s.U[ider][gr.iU[i]] -= y[gr.jU[i]] * gr.scaleU[i]; end end
    if ider==1          for i ∈ eachindex(gr.iA); s.A[      gr.iA[i]] -= y[gr.jA[i]] * gr.scaleA[i]; end end
end
function increment!(s::State,ider::𝕫,y::AbstractVector{𝕣},gr::DofGroup) 
    if ider≤length(s.Λ) for i ∈ eachindex(gr.iΛ); s.Λ[ider][gr.iΛ[i]] += y[gr.jΛ[i]] * gr.scaleΛ[i]; end end
    if ider≤length(s.X) for i ∈ eachindex(gr.iX); s.X[ider][gr.iX[i]] += y[gr.jX[i]] * gr.scaleX[i]; end end
    if ider≤length(s.U) for i ∈ eachindex(gr.iU); s.U[ider][gr.iU[i]] += y[gr.jU[i]] * gr.scaleU[i]; end end
    if ider==1          for i ∈ eachindex(gr.iA); s.A[      gr.iA[i]] += y[gr.jA[i]] * gr.scaleA[i]; end end
end
function set!(s::State,ider::𝕫,y::AbstractVector{𝕣},gr::DofGroup) 
    if ider≤length(s.Λ) s.Λ[ider] .= 0 end
    if ider≤length(s.X) s.X[ider] .= 0 end
    if ider≤length(s.U) s.U[ider] .= 0 end
    if ider==1          s.A       .= 0 end
    if ider≤length(s.Λ) for i ∈ eachindex(gr.iΛ); s.Λ[ider][gr.iΛ[i]] = y[gr.jΛ[i]] * gr.scaleΛ[i]; end end
    if ider≤length(s.X) for i ∈ eachindex(gr.iX); s.X[ider][gr.iX[i]] = y[gr.jX[i]] * gr.scaleX[i]; end end
    if ider≤length(s.U) for i ∈ eachindex(gr.iU); s.U[ider][gr.iU[i]] = y[gr.jU[i]] * gr.scaleU[i]; end end
    if ider==1          for i ∈ eachindex(gr.iA); s.A[      gr.iA[i]] = y[gr.jA[i]] * gr.scaleA[i]; end end
end
function get!(s::State,ider::𝕫,y::AbstractVector{𝕣},gr::DofGroup) 
    if ider≤length(s.Λ)for i ∈ eachindex(gr.iΛ); y[gr.jΛ[i]] = s.Λ[ider][gr.iΛ[i]] / gr.scaleΛ[i]; end end
    if ider≤length(s.X)for i ∈ eachindex(gr.iX); y[gr.jX[i]] = s.X[ider][gr.iX[i]] / gr.scaleX[i]; end end
    if ider≤length(s.U)for i ∈ eachindex(gr.iU); y[gr.jU[i]] = s.U[ider][gr.iU[i]] / gr.scaleU[i]; end end
    if ider==1         for i ∈ eachindex(gr.iA); y[gr.jA[i]] = s.A[      gr.iA[i]] / gr.scaleA[i]; end end
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
nodofs(    model::Model,dis) = DofGroup(dis, 𝕫[],𝕫[],𝕫[],𝕫[])
allΛdofs(  model::Model,dis) = DofGroup(dis, 1:getndof(model,:X),𝕫[],𝕫[],𝕫[])
allXdofs(  model::Model,dis) = DofGroup(dis, 𝕫[],1:getndof(model,:X),𝕫[],𝕫[])
allUdofs(  model::Model,dis) = DofGroup(dis, 𝕫[],𝕫[],1:getndof(model,:U),𝕫[])
allAdofs(  model::Model,dis) = DofGroup(dis, 𝕫[],𝕫[],𝕫[],1:getndof(model,:A))
allXUdofs( model::Model,dis) = DofGroup(dis, 𝕫[],1:getndof(model,:X),1:getndof(model,:U),𝕫[])
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

# asm[iarray,ieletyp][ieledof/i,iele] has value zero for terms from element gradient/hessian that are not to be added in. Otherwise, the value they
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
function gradientpartition(nΛ,nX,nU,nA,nXder=1,nUder=nXder)
    iΛ          =                           (1:nΛ)
    iX          = nΛ                      .+(1:nX)
    iU          = nΛ           +(nX*nXder).+(1:nU) 
    iA          = nΛ+(nX*nXder)+(nU*nUder).+(1:nA)
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
                i = ieledof+nieledof*(jeledof-1)
                asm[ieletyp][i,iele] = imoddof+nimoddof*(jmoddof-1)
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
                i = ieledof+nieledof*(jeledof-1) 
                asm[ieletyp][i,iele] = K[ipair]  
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
# assemble_!         Muscade/Assemble.jl      loop over elements within type (typestable) (except Acost: do nothing)
# assemble!          Muscade/Assemble.jl      loop over element types (barrier function) 
# assemble_!         Muscade/Assemble.jl      loop over elements within type (typestable) if Acost otherwise do nothing if Acost
# addin!             Muscade/SomeSolver.jl    do adiff and add-in, solver-specific
# getresidual        Muscade/Assemble.jl      typechecking the call from addin! call residual or Lagrangian as available, check for NaNs
# residual           MyElement.jl             the element code 


abstract type Assembly end # solver define concrete "assemblies" which is a collection of matrices and solvers wanted for a phase in the solution process
struct assemble!{mission} end
struct assemble_!{mission} end
struct assembleA!{mission} end
struct assembleA_!{mission} end
struct addin!{mission} end

# standard assembler, everything except time-independant A-cost
# Assemble a state using the values of variables in the state.
function assemble!{mission}(out::Assembly,asm,dis,model,state,Δt,dbg) where{mission}
    zero!(out)
    for ieletyp = 1:lastindex(model.eleobj)
        eleobj  = model.eleobj[ieletyp]
        assemble_!{mission}(out,view(asm,:,ieletyp),dis.dis[ieletyp],eleobj,state,state.time,Δt,state.SP,(dbg...,ieletyp=ieletyp))
    end
end
assemble_!{         mission}(out::Assembly,asm,dis,eleobj::Vector{<:Acost},state::State{nΛder,nXder,nUder},t,Δt,SP,dbg) where{mission,nΛder,nXder,nUder} = nothing
function assemble_!{mission}(out::Assembly,asm,dis,eleobj                 ,state::State{nΛder,nXder,nUder},t,Δt,SP,dbg) where{mission,nΛder,nXder,nUder}
    for iele  = 1:lastindex(eleobj)
        index = dis.index[iele]
        Λe    = NTuple{nΛder}(λ[index.X] for λ∈state.Λ)
        Xe    = NTuple{nXder}(x[index.X] for x∈state.X)
        Ue    = NTuple{nUder}(u[index.U] for u∈state.U)
        Ae    =         state.A[index.A]
        addin!{mission}(out,asm,iele,dis.scale,eleobj[iele],Λe,Xe,Ue,Ae, t,Δt,SP,(dbg...,iele=iele)) # defined by solver.  Called for each element. But the asm that is passed
    end                                                                              # is of the form asm[iarray][i,iele], because addin! will add to all arrays in one pass
end

# TODO stateless assembler.  Who uses this??? 
# if we keep this: make standard assemble! above call the assemble_! below
#
# function assemble!{mission}(out::Assembly,asm,dis,model,Λ,X,U,A,t,Δt,dbg) where{mission}
#     zero!(out)
#     for ieletyp = 1:lastindex(model.eleobj)
#         eleobj  = model.eleobj[ieletyp]
#         assemble_!{mission}(out,view(asm,:,ieletyp),dis.dis[ieletyp],eleobj,Λ,X,U,A,t,Δt,nothing,(dbg...,ieletyp=ieletyp))
#     end
# end
# assemble_!{mission}(out::Assembly,asm,dis,eleobj::Acost,Λ,X,U,A,t,Δt,SP,dbg) where{mission} = nothing
# function assemble_!{mission}(out::Assembly,asm,dis,eleobj,Λ::NTuple{nΛder},X::NTuple{nXder},U::NTuple{nUder},A,t,Δt,SP,dbg) where{mission,nΛder,nXder,nUder}
#     for iele  = 1:lastindex(eleobj)
#         index = dis.index[iele]
#         Λe    = NTuple{nΛder}(λ[index.X] for λ∈Λ)
#         Xe    = NTuple{nXder}(x[index.X] for x∈X)
#         Ue    = NTuple{nUder}(u[index.U] for u∈U)
#         Ae    = A[index.A]
#         addin!{mission}(out,asm,iele,dis.scale,eleobj[iele],Λe,Xe,Ue,Ae, t,Δt,SP,(dbg...,iele=iele)) # defined by solver.  Called for each element. But the asm that is passed
#     end                                                                              # is of the form asm[iarray][i,iele], because addin! will add to all arrays in one pass
# end

# assemble the time-independant A-cost
function assembleA!{mission}(out::Assembly,asm,dis,model,state,dbg) where{mission}
    zero!(out)
    for ieletyp = 1:lastindex(model.eleobj)
        eleobj  = model.eleobj[ieletyp]
        assembleA_!{mission}(out,view(asm,:,ieletyp),dis.dis[ieletyp],eleobj,state,(dbg...,assembleA=true,ieletyp=ieletyp))
    end
end
assembleA_!{mission}(out::Assembly,asm,dis,eleobj,state,dbg) where{mission} = nothing
function assembleA_!{mission}(out::Assembly,asm,dis,eleobj::Vector{A},state,dbg) where{mission,A<:Acost}
    for iele  = 1:lastindex(eleobj)
        Ae    = state.A[dis.index[iele].A]
        addin!{mission}(out,asm,iele,dis.scale,eleobj[iele],Ae,(dbg...,iele=iele))                
    end
end


############# Tools for addin!



#### zero!
"""
    zero!(a)

Set to zero all elements of an arrays. If `a` is sparse, 
the vector `nzval` of values is set to zero and the sparsity structure is unchanged.
""" 
function zero!(out::Base.RefValue)
    out[] = 0
end
function zero!(out::AbstractArray)
    for i∈eachindex(out)
        out[i] = 0
    end
end
function zero!(out::AbstractSparseArray)
    for i∈eachindex(out.nzval)
        out.nzval[i] = 0
    end
end
function zero!(out::Tuple)
    for oᵢ ∈ out
        zero!(oᵢ)
    end
end

#### extract value or derivatives from a SVector 'a' of adiffs, and add it directly into vector, full matrix or nzval of sparse matrix 'out'.
# given a vector a
# ia:    where in the vector to pick data
# ida:   where in the partials to pick data
# iasm:  where in first dim of out to put data      (if transpose: second)
# idasm: where in the second dim of out to put data (if transpose: first)
#
# out[asm[:   ,iele]] += a
# out[asm[iasm,iele]] += a      # build: 'a' is only a part of the element vector (FreqXU)   
# out[asm[:,   iele]] += a[ia]  # split: parts of 'a' are assembled (DirectXUA)   
# out[asm[iasm,iele]] += a[ia]  # not used
function add_value!(out::𝕣1,asm,iele,a::SVector{Na,<:ℝ},ia=1:Na ; iasm=idvec,Δt=idmult) where{Na}
    for (i,iaᵢ) ∈ enumerate(ia)
        iout = asm[iasm[i],iele]
        if iout≠0 
            out[iout]+=VALUE(a[iaᵢ])*Δt 
        end
    end
end   

struct   add_∂!{P,S,T} end # to allow syntax with type-parameter P: precedence, S: :plus|:minus, T: :transpose|:notranspose
#### derivatives from a SVector 'a' of adiffs, and add it directly into full matrix or nzval of sparse matrix 'out'.
# given a vector a
# ia:    where in the vector to pick data
# ida:   where in the partials to pick data
# iasm:  where in first dim of out to put data      (if transpose: second)
# idasm: where in the second dim of out to put data (if transpose: first)
#
# ∀ i,j out[asm[iasm[i]+length(ia)*idasm[j]-1,iele]] += a[ia[i]].dx[ida[j]]  
#
# out::Array is deliberately vague, encompassing 𝕣2 and sparse.nzval
function add_∂!{P,S,T}(out::Array,asm,iele,a::SVector{Na,∂ℝ{P,Nda,R}},ia=1:Na,ida=1:Nda ; iasm=idvec,idasm=idvec,Δt=idmult) where{P,Nda,R,Na,S,T}
    for (i,iaᵢ) ∈ enumerate(ia), (j,idaⱼ) ∈ enumerate(ida)
        k = if T==:transpose   idasm[j]+length(ida)*( iasm[i]-1)   
        elseif T==:notranspose iasm[ i]+length( ia)*(idasm[j]-1)  
        else   muscadeerror((;T=T),"Illegal value of parameter T")    
        end
        iout = asm[k,iele]
        if iout≠0
            if     S==:plus   out[iout]+=a[iaᵢ].dx[idaⱼ]*Δt  
            elseif S==:minus  out[iout]-=a[iaᵢ].dx[idaⱼ]*Δt  
            else   muscadeerror((;S=S),"Illegal value of parameter S")    
            end
        end
    end
end  
add_∂!{P    }(                                     args...;kwargs...) where{P         } = add_∂!{P,:plus             }(args...;kwargs...) 
add_∂!{P,S  }(                                     args...;kwargs...) where{P,S       } = add_∂!{P,S    ,:notranspose}(args...;kwargs...) 
add_∂!{P,S,T}(out::SparseMatrixCSC,                args...;kwargs...) where{P,S,T     } = add_∂!{P,S    ,T           }(out.nzval, args...;kwargs...)
add_∂!{P,S,T}(out::Array,asm,iele,a::SVector{Na,R},args...;kwargs...) where{P,S,T,Na,R} = nothing # if P does not match

function add_∂!{P,S}(out::Vector,asm, iele, a::SVector{Na,∂ℝ{P,Nda,R}},ia,ida::𝕫 ; Δt=idmult) where{P,S,Nda,R,Na} # Lλr::Vector in Newmark-β context
    for (i,iaᵢ) ∈ enumerate(ia)
        iout = asm[i,iele]
        if iout≠0
            if     S==:plus   out[iout]+=a[iaᵢ].dx[ida]*Δt  
            elseif S==:minus  out[iout]-=a[iaᵢ].dx[ida]*Δt  
            else   muscadeerror((;S=S),"Illegal value of parameter S")    
            end
        end
    end
end   
function add_∂!{P,S,T}(out::𝕣0,a::SVector{Na,∂ℝ{P,Nda,R}},ia::𝕫,ida::𝕫;Δt=idmult) where{P,S,T,Nda,R,Na} # Lrr, scalar in Newmark-β context
    if     S==:plus   out[]+=a[ia].dx[ida]*Δt  
    elseif S==:minus  out[]-=a[ia].dx[ida]*Δt  
    else   muscadeerror((;S=S),"Illegal value of parameter S")    
    end
end

####### called by addin!, and by nested elements to "get a Lagrangian" and "get a residual"
# 1) comprehensive check of the types of arguments, to help catch bugs in solvers and elements at compile time
# 2) if a residual is wanted by the solver and only Lagrangian is implemented by the element (or the other way around), do some magic
# 3) check for NaNs in the results 
#
# Note that getLagrangian receives Λ::SVector. addin! by contrast receives Λ::NTuple{SVector}, this is not a bug


## Type unstable(?) stopgap.  Solution? Get @espy to generate hasresidual(  ::Eleobj)=true, with false fallback
hasresidual(  ::Eleobj) where{Eleobj} = Val(hasmethod(residual  ,(Eleobj,       NTuple,NTuple,𝕣1,𝕣,NamedTuple,NamedTuple)))
haslagrangian(::Eleobj) where{Eleobj} = Val(hasmethod(lagrangian,(Eleobj,NTuple,NTuple,NTuple,𝕣1,𝕣,NamedTuple,NamedTuple)))

# dispatcher
function getresidual(eleobj::Eleobj,X,U,A,t,SP,dbg,req) where{Eleobj} 
    R,FB,eleres = getresidual(eleobj,hasresidual(eleobj),haslagrangian(eleobj),no_second_order(Eleobj),X,U,A,t,SP,dbg,req) 
    R isa SVector{length(X[1])} || muscadeerror(dbg,@sprintf("R, from residual(%s,...), must be a SVector of length nXdof",Eleobj))
    hasnan(R,FB) && muscadeerror((dbg...,t=t,SP=SP ),@sprintf("residual(%s,...) returned NaN in R, FB or derivatives",Eleobj))  
    return R,FB,eleres
end
function getresidual(eleobj::Eleobj,X,U,A,t,SP,dbg,req...) where{Eleobj} 
    R,FB     = getresidual(eleobj,hasresidual(eleobj),haslagrangian(eleobj),no_second_order(Eleobj),X,U,A,t,SP,dbg    ) 
    R isa SVector{length(X[1])} || muscadeerror(dbg,@sprintf("R, from residual(%s,...), must be a SVector of length nXdof",Eleobj))
    hasnan(R,FB) && muscadeerror((dbg...,t=t,SP=SP ),@sprintf("residual(%s,...) returned NaN in R, FB or derivatives",Eleobj))  
    return R,FB
end

# has residual
function getresidual(eleobj::Eleobj,hasres::Val{true},haslag,nso::Val{false}, X::NTuple{Ndx,SVector{Nx}}, 
        U::NTuple{Ndu,SVector{Nu}}, A::SVector{Na} ,t::ℝ,SP,dbg,req)     where{Eleobj<:AbstractElement,Ndx,Nx,Ndu,Nu,Na} 
    R,FB,eleres = residual(  eleobj,  X,U,A,t,SP,dbg,req)
end
function getresidual(eleobj::Eleobj,hasres::Val{true},haslag,nso::Val{false}, X::NTuple{Ndx,SVector{Nx}}, 
        U::NTuple{Ndu,SVector{Nu}}, A::SVector{Na} ,t::ℝ,SP,dbg    )     where{Eleobj<:AbstractElement,Ndx,Nx,Ndu,Nu,Na} 
    R,FB       = residual(  eleobj,  X,U,A,t,SP,dbg    )
end

# has residual with no second order
function getresidual(eleobj::Eleobj,hasres::Val{true},haslag,nso::Val{true}, X::NTuple{Ndx,SVector{Nx}}, 
        U::NTuple{Ndu,SVector{Nu}}, A::SVector{Na} ,t::ℝ,SP,dbg,req)     where{Eleobj<:AbstractElement,Ndx,Nx,Ndu,Nu,Na} 
    P,N = precedence((X,U,A,t)),npartial((X,U,A,t)) 
    X1,U1,A1,t1 = to_order{min(1,P),N}((X,U,A,t))
    R,FB,eleres = residual(  eleobj,  X1,U1,A1,t1,SP,dbg,req)
    return to_order{P,N}(R),FB,to_order{P,N}(eleres) 
end
function getresidual(eleobj::Eleobj,hasres::Val{true},haslag,nso::Val{true}, X::NTuple{Ndx,SVector{Nx}}, 
        U::NTuple{Ndu,SVector{Nu}}, A::SVector{Na} ,t::ℝ,SP,dbg    )     where{Eleobj<:AbstractElement,Ndx,Nx,Ndu,Nu,Na} 
    P,N = precedence((X,U,A,t)),npartial((X,U,A,t)) 
    X1,U1,A1,t1 = to_order{min(1,P),N}((X,U,A,t))
    R,FB        = residual(  eleobj,  X1,U1,A1,t1,SP,dbg    )
    return to_order{P,N}(R),FB 
end

# has lagrangian
function getresidual(eleobj::Eleobj,hasres::Val{false},haslag::Val{true},nso, X::NTuple{Ndx,SVector{Nx}}, 
        U::NTuple{Ndu,SVector{Nu}}, A::SVector{Na} ,t::ℝ,SP,dbg,req)     where{Eleobj<:AbstractElement,Ndx,Nx,Ndu,Nu,Na} 
    λ            = SVector{Nx,𝕣}(0. for i=1:Nx)
    Pλ,_,Λ       = variate(λ,context=(X,U,A,t))
    L,FB,eleres  = lagrangian(eleobj,Λ,X,U,A,t,SP,dbg,req)    
    R            = ∂{Pλ,Nx}(L)
    return R,FB,eleres
end
function getresidual(eleobj::Eleobj,hasres::Val{false},haslag::Val{true},nso, X::NTuple{Ndx,SVector{Nx}}, 
        U::NTuple{Ndu,SVector{Nu}}, A::SVector{Na} ,t::ℝ,SP,dbg   )     where{Eleobj<:AbstractElement,Ndx,Nx,Ndu,Nu,Na} 
    λ            = SVector{Nx,𝕣}(0. for i=1:Nx)
    Pλ,_,Λ       = variate(λ,context=(X,U,A,t))
    L,FB         = lagrangian(eleobj,Λ,X,U,A,t,SP,dbg   )    
    R            = ∂{Pλ,Nx}(L)
    return R,FB
end

# has nothing !?!
getresidual(eleobj::Eleobj,hasres::Val{false},haslag::Val{false},nso, X,U,A,t,SP,dbg,req...) where{Eleobj} =
    muscadeerror((dbg...,t=t,SP=SP),@sprintf("Element %s must have method 'Muscade.lagrangian' or/and 'Muscade.residual' with correct interface",Eleobj))

# Acost
function getlagrangian(eleobj::Acost, A::SVector, SP, dbg,req) 
    L,FB,eleres = lagrangian(eleobj,A,SP,dbg,req)    
    L isa ℝ || muscadeerror(dbg,@sprintf("L, from lagrangian(%s,...), must be a ℝ, alias Real",Eleobj))
    hasnan(L,FB) && muscadeerror((dbg...,t=t,SP=SP),"lagrangian(Acost,...) returned NaN in L, FB or derivatives") 
    return L,FB,eleres  
end
function getlagrangian(eleobj::Acost, A::SVector, SP, dbg) 
    L,FB = lagrangian(eleobj,A,SP,dbg)    
    L isa ℝ || muscadeerror(dbg,@sprintf("L, from lagrangian(%s,...), must be a ℝ, alias Real",Eleobj))
    hasnan(L,FB) && muscadeerror((dbg...,t=t,SP=SP),@sprintf("lagrangian(Acost,...) returned NaN in L, FB or derivatives",Eleobj)) 
    return L,FB  
end

# dispatcher
function getlagrangian(eleobj::Eleobj, Λ::SVector{Nx}, X::NTuple{Ndx,SVector{Nx}}, 
        U::NTuple{Ndu,SVector{Nu}}, A::SVector{Na}, t::ℝ,SP,dbg,req)     where{Eleobj<:AbstractElement,Ndx,Nx,Ndu,Nu,Na} 
    L,FB,eleres = getlagrangian(eleobj,hasresidual(eleobj),haslagrangian(eleobj),no_second_order(Eleobj),Λ,X,U,A,t,SP,dbg,req)    
    L isa ℝ || muscadeerror(dbg,@sprintf("L, from lagrangian(%s,...), must be a ℝ, alias Real",Eleobj))
    hasnan(L,FB) && muscadeerror((dbg...,t=t,SP=SP),@sprintf("lagrangian(%s,...) returned NaN in L, FB or derivatives",Eleobj)) 
    return L,FB,eleres  
end
function getlagrangian(eleobj::Eleobj, Λ::SVector{Nx}, X::NTuple{Ndx,SVector{Nx}}, 
        U::NTuple{Ndu,SVector{Nu}}, A::SVector{Na}, t::ℝ,SP,dbg,req...)     where{Eleobj<:AbstractElement,Ndx,Nx,Ndu,Nu,Na} 
    L,FB           = getlagrangian(eleobj,hasresidual(eleobj),haslagrangian(eleobj),no_second_order(Eleobj),Λ,X,U,A,t,SP,dbg       )    
    L isa ℝ || muscadeerror(dbg,@sprintf("L, from lagrangian(%s,...), must be a ℝ, alias Real",Eleobj))
    hasnan(L,FB) && muscadeerror((dbg...,t=t,SP=SP),@sprintf("lagrangian(%s,...) returned NaN in L, FB or derivatives",Eleobj)) 
    return L,FB        
end

# has lagrangian
function getlagrangian(eleobj::Eleobj,hasres,haslag::Val{true},nso, Λ::SVector{Nx}, X::NTuple{Ndx,SVector{Nx}}, 
        U::NTuple{Ndu,SVector{Nu}}, A::SVector{Na}, t::ℝ,SP,dbg,req)     where{Eleobj<:AbstractElement,Ndx,Nx,Ndu,Nu,Na} 
    L,FB,eleres = lagrangian(eleobj,Λ,X,U,A,t,SP,dbg,req)
end
function getlagrangian(eleobj::Eleobj,hasres,haslag::Val{true},nso, Λ::SVector{Nx}, X::NTuple{Ndx,SVector{Nx}}, 
        U::NTuple{Ndu,SVector{Nu}}, A::SVector{Na}, t::ℝ,SP,dbg   )     where{Eleobj<:AbstractElement,Ndx,Nx,Ndu,Nu,Na} 
    L,FB        = lagrangian(eleobj,Λ,X,U,A,t,SP,dbg    )
end

# has residual
function getlagrangian(eleobj::Eleobj,hasres::Val{true},haslag::Val{false},nso::Val{false}, Λ::SVector{Nx}, X::NTuple{Ndx,SVector{Nx}}, 
        U::NTuple{Ndu,SVector{Nu}}, A::SVector{Na}, t::ℝ,SP,dbg,req)     where{Eleobj<:AbstractElement,Ndx,Nx,Ndu,Nu,Na} 
    R,FB,eleres = residual(  eleobj,  X,U,A,t,SP,dbg,req)
    L = Λ ∘₁ R
    return L,FB,eleres 
end
function getlagrangian(eleobj::Eleobj,hasres::Val{true},haslag::Val{false},nso::Val{false}, Λ::SVector{Nx}, X::NTuple{Ndx,SVector{Nx}}, 
        U::NTuple{Ndu,SVector{Nu}}, A::SVector{Na}, t::ℝ,SP,dbg    )     where{Eleobj<:AbstractElement,Ndx,Nx,Ndu,Nu,Na} 
    R,FB           = residual(  eleobj,  X,U,A,t,SP,dbg    )
    L = Λ ∘₁ R
    return L,FB           
end

# has residual with no second order
function getlagrangian(eleobj::Eleobj,hasres::Val{true},haslag::Val{false},nso::Val{true}, Λ::SVector{Nx}, X::NTuple{Ndx,SVector{Nx}}, 
        U::NTuple{Ndu,SVector{Nu}}, A::SVector{Na}, t::ℝ,SP,dbg,req)     where{Eleobj<:AbstractElement,Ndx,Nx,Ndu,Nu,Na} 
    P,N = precedence((X,U,A,t)),npartial((X,U,A,t)) 
    X1,U1,A1,t1 = to_order{min(1,P),N}((X,U,A,t)) 
    R,FB,eleres = residual(  eleobj,  X1,U1,A1,t1,SP,dbg,req)
    L           = Λ ∘₁ to_order{P,N}(R) # to avoid loosing symmetry of Hessian...
    return L,FB,to_order{P,N}(eleres) 
end
function getlagrangian(eleobj::Eleobj,hasres::Val{true},haslag::Val{false},nso::Val{true}, Λ::SVector{Nx}, X::NTuple{Ndx,SVector{Nx}}, 
        U::NTuple{Ndu,SVector{Nu}}, A::SVector{Na}, t::ℝ,SP,dbg    )     where{Eleobj<:AbstractElement,Ndx,Nx,Ndu,Nu,Na} 
    P,N = precedence((X,U,A,t)),npartial((X,U,A,t)) 
    X1,U1,A1,t1 = to_order{min(1,P),N}((X,U,A,t)) 
    R,FB        = residual(  eleobj,  X1,U1,A1,t1,SP,dbg    )
    L           = Λ ∘₁ to_order{P,N}(R) # to avoid loosing symmetry of Hessian...
    return L,FB 
end

# has nothing !?!
getlagrangian(eleobj::Eleobj,hasres::Val{false},haslag::Val{false},nso, Λ,X,U,A,t,SP,dbg,req...)     where{Eleobj} =
    muscadeerror((dbg...,t=t,SP=SP),@sprintf("Element %s must have method 'Muscade.lagrangian' or/and 'Muscade.residual' with correct interface",Eleobj))
