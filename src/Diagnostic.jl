
##################  describe
"""
    describe(model,spec)

Print out information about `model`.
`spec` can be 
- an `EleID` to describe an element,
- a `DofID` to describe a dof.
- a `NodID` to describe a node,
- `:doftyp` to obtain a list of doftypes, 
- `:dof` to obtain a list of dofs or 
- `:eletyp` for a list of element types.

    describe(state,[class=:all])

Provide a description of the dofs stored in `state`.
`class` can be either `:all`, `:Λ`, `:ΛX`, `:X`, `:U`, `:A` or `:scale`.

See also: [`addelement!`](@ref), [`addnode!`](@ref)
"""
function describe(model::Model,eleID::EleID)
    try 
        dof = model.ele[eleID] 
    catch
        printstyled("Not a valid EleID\n",color=:red,bold=true)
        return
    end
    e  = model.ele[eleID]
    eo = model.eleobj[eleID]
    @printf "Element with EleID(%i,%i)\n" eleID.ieletyp eleID.iele 
    @printf "   model.eleobj[%i][%i]::" eleID.ieletyp eleID.iele 
    printstyled(@sprintf("%s\n",typeof(eo)),color=:cyan)
    @printf "   model.ele[%i][%i]:\n" eleID.ieletyp eleID.iele
    for dofid ∈ e.dofID
        dof    = model.dof[dofid]
        nod    = model.nod[dof.nodID]
        doftyp = model.doftyp[dof.idoftyp]
        @printf "      NodID(%i), class=:%s, field=:%-12s\n" dof.nodID.inod doftyp.class doftyp.field 
    end
end
function describe(model::Model,dofID::DofID)
    try 
        dof = model.dof[dofID] 
    catch
        printstyled("Not a valid DofID\n",color=:red,bold=true)
        if dofID.class==:Λ
            @printf "Optimisation solvers introduce a one-to-one correspondance between :Λ-dofs and :X-dofs, \nbut :Λ-dofs are not part of the model description: try DofID(:X,...)\n"
        end
        return
    end
    dof     = model.dof[dofID] 
    doftyp  = model.doftyp[dof.idoftyp]
    @printf "Degree of freedom with DofID(:%s,%i)\n" dofID.class dofID.idof
    @printf "   model.dof.%s[%i]:\n" dofID.class dofID.idof
    @printf "   NodID(%i), class=:%s, field=:%-12s\n" dof.nodID.inod dofID.class doftyp.field 
    @printf "   elements:\n"
    for eleid ∈ dof.eleID
        @printf "      EleID(%i,%i), " eleid.ieletyp eleid.iele 
        printstyled(@sprintf("%s\n",eltype(model.eleobj[eleid.ieletyp])),color=:cyan)
    end
    if dofID.class == :X
        @printf "   Output in state[istep].X[ider+1][%i] and state[istep].Λ[%i]\n" dofID.idof dofID.idof    
    elseif dofID.class ==:U
            @printf "   Output in state[istep].U[ider][%i]\n" dofID.idof   
        elseif dofID.class == :A
        @printf "   Output in state[istep].A[%i]\n" dofID.idof   
    end            
end
function describe(model::Model,nodID::NodID)
    try 
        nod = model.nod[nodID] 
    catch
        printstyled("Not a valid NodID\n",color=:red,bold=true)
        return
    end
    nod = model.nod[nodID]
    @printf "Node with NodID(%i)\n" nodID.inod
    @printf "   model.nod[%i]:\n" nodID.inod
    nc = length(nod.coord)
    @printf "   coord=[" 
    for ic=1:nc-1
        @printf "%g," nod.coord[ic] 
    end
    if nc>0
        @printf "%g" nod.coord[nc] 
    end
    @printf "]\n" 
    @printf "   dof (degrees of freedom):\n"
    for dofid ∈ nod.dofID
        dof = model.dof[dofid]
        doftyp = model.doftyp[dof.idoftyp]
        @printf "      DofID(:%s,%i), class=:%s, inod=%i, field=:%-12s\n" dofid.class dofid.idof dofid.class nodID.inod doftyp.field    
    end
    @printf "   elements:\n"
    for eleID ∈ nod.eleID
        @printf "      EleID(%i,%i), " eleID.ieletyp eleID.iele 
        printstyled(@sprintf("%s\n",typeof(model.eleobj[eleID])),color=:cyan)
    end
 end
function describe(model::Model,s::Symbol)
    @printf "\nModel '%s'\n" model.ID
    if s==:dof
        for class ∈ (:X,:U,:A)
            ndof    = getndof(model,class)
            ndof>0 && @printf "\n   Dofs of class :%s\n" class
            for idof = 1:ndof
                dof     = model.dof[class][idof] 
                doftyp  = model.doftyp[dof.idoftyp]
                @printf "      %6d. field= :%-15s NodID(%i)\n" idof doftyp.field dof.nodID.inod 
            end
        end
    elseif s==:doftyp
        for doftyp ∈ model.doftyp
            doftyp.class==:X && @printf "      class= :%-5s field= :%-15s\n" doftyp.class doftyp.field 
        end
        for doftyp ∈ model.doftyp
            doftyp.class==:U && @printf "      class= :%-5s field= :%-15s\n" doftyp.class doftyp.field 
        end
        for doftyp ∈ model.doftyp
            doftyp.class==:A && @printf "      class= :%-5s field= :%-15s\n" doftyp.class doftyp.field 
        end
    elseif s==:eletyp
        et = eletyp(model)
        for i∈eachindex(et)
            @printf "    %6d elements of type %s\n" length(model.eleobj[i]) et[i]
        end
    end
end
function describeX(state::State)
    model = state.model
    nX    = getndof(model,:X)
    nder  = length(state.X)
    for iX = 1:nX
        dofID   = DofID(:X,iX)
        dof     = model.dof[dofID] 
        @printf "NodID(%i), class=:%s, field=:%-15s   " dof.nodID.inod dofID.class state.dis.fieldX[iX]
        for ider = 1:nder
            @printf "%15g " state.X[ider][iX]
        end
        @printf "\n" 
    end
end
function describeΛX(state::State)
    model = state.model
    nX    = getndof(model,:X)
    nXder = length(state.X)
    nΛder = length(state.Λ)
    for iX = 1:nX
        dofID   = DofID(:X,iX)
        dof     = model.dof[dofID] 
        @printf "NodID(%i), class=:%s, field=:%-15s, " dof.nodID.inod dofID.class state.dis.fieldX[iX] 
        @printf "X="
        for ider = 1:nXder
            @printf "%15g " state.X[ider][iX]
        end
        @printf ", Λ="
        for ider = 1:nΛder
            @printf "%15g " state.Λ[ider][iX]
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
        @printf "NodID(%i), class=:%s, field=:%-15s   " dof.nodID.inod dofID.class state.dis.fieldU[iU]
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
        @printf "NodID(%i), class=:%s, field=:%-15s   %15g\n" dof.nodID.inod dofID.class state.dis.fieldA[iA] state.A[iA] 
    end
end
function describeScale(state::State)
    dis = state.dis
    bigΛ            = Dict{Symbol,𝕣}()
    for idof        ∈ eachindex(state.Λ)
        field       = dis.fieldX[idof]
        bigΛ[field] = max(get(bigΛ,field,0.),abs(state.Λ[1][idof]))
    end
    for field       ∈ keys(bigΛ)
        @printf "class= :Λ field= :%-15s  max(|dof|)= %g\n" field bigΛ[field]
    end
    bigX            = Dict{Symbol,𝕣}()
    for idof        ∈ eachindex(state.X[1])
        field       = dis.fieldX[idof]
        bigX[field] = max(get(bigX,field,0.),abs(state.X[1][idof]))
    end
    for field       ∈ keys(bigX)
        @printf "class= :X field= :%-15s  max(|dof|)= %g\n" field bigX[field]
    end
    bigU            = Dict{Symbol,𝕣}()
    for idof        ∈ eachindex(state.U[1])
        field       = dis.fieldU[idof]
        bigU[field] = max(get(bigU,field,0.),abs(state.U[1][idof]))
    end
    for field       ∈ keys(bigU)
        @printf "class= :U field= :%-15s  max(|dof|)= %g\n" field bigU[field]
    end
    bigA            = Dict{Symbol,𝕣}()
    for idof        ∈ eachindex(state.A)
        field       = dis.fieldA[idof]
        bigA[field] = max(get(bigA,field,0.),abs(state.A[idof]))
    end
    for field       ∈ keys(bigA)
        @printf "class= :A field= :%-15s  max(|dof|)= %g\n" field bigA[field]
    end
end 

function getfield_(dg::DofGroup,class)
    return if class == :Λ dg.fieldΛ
    elseif    class == :X dg.fieldX
    elseif    class == :U dg.fieldU
    elseif    class == :A dg.fieldA
    else muscadeerror("Class must be :Λ, :X,:U or :A")
    end 
end    
function getj(dg::DofGroup,class)
    return if class == :Λ dg.jΛ
    elseif    class == :X dg.jX
    elseif    class == :U dg.jU
    elseif    class == :A dg.jA
    else muscadeerror("Class must be :Λ, :X,:U or :A")
    end 
end    
function describeScale(v::Vector,dofgr::DofGroup) # Actually for debugging use
    big             = Dict{Tuple{Symbol, Symbol},𝕣}()
    for class ∈ (:Λ,:X,:U,:A)
        js              = getj(    dofgr,class)
        fields          = getfield_(dofgr,class)
        for (i,field)   ∈ enumerate(fields)
            big[(class,field)] = max(get(big,(class,field),0.),abs(v[js[i]]))
        end
    end
    for (class,field) ∈ keys(big)
        @printf "|v[%s-%s]| ≤ %g\n" class field big[(class,field)]
    end
end
function describeScale(m::AbstractMatrix,idofgr::DofGroup,jdofgr::DofGroup=idofgr) # Actually for debugging use
    big             = Dict{Tuple{Symbol, Symbol,Symbol, Symbol},𝕣}()
    for iclass ∈ (:Λ,:X,:U,:A)
        ijs              = getj(    idofgr,iclass)
        ifields          = getfield_(idofgr,iclass)
        for jclass ∈ (:Λ,:X,:U,:A)
            jjs              = getj(    jdofgr,jclass)
            jfields          = getfield_(jdofgr,jclass)
            for (i,ifield)   ∈ enumerate(ifields)
                for (j,jfield) ∈ enumerate(jfields)
                    key           = (iclass,ifield,jclass,jfield)
                    val           = abs(m[ijs[i],jjs[j]])
                    if val>0
                        big[key]      = max(get(big,key,0.),abs(m[ijs[i],jjs[j]]))
                    end
                end
            end
        end
    end
    for (iclass,ifield,jclass,jfield) ∈ keys(big)
        @printf "|m[%s-%s,%s-%s]| ≤ %g\n" iclass ifield jclass jfield big[(iclass,ifield,jclass,jfield)]
    end
end

function describe(state::State{nΛder,nXder,nUder,TSP},class::Symbol=:all) where{nΛder,nXder,nUder,TSP}
    @printf("State{nΛder=%i,nXder=%i,nUder=%i,TSP=%s}\n",nΛder,nXder,nUder,TSP)
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
    elseif class == :scale 
        describeScale(state)    
    else
        printstyled("Not a valid class\n",color=:red,bold=true)
    end
end




#################### study_scale


mutable struct Assemblystudy_scale{Tz,Tzz}  <:Assembly
    Lz    :: Tz
    Lzz   :: Tzz
end   
function prepare(::Type{Assemblystudy_scale},model,dis) 
    dofgr              = allΛXUAdofs(model,dis)
    nZ                 = getndof(dofgr)
    narray,neletyp     = 2,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    Lz                 = asmvec!(view(asm,1,:),dofgr,dis) 
    Lzz                = asmmat!(view(asm,2,:),view(asm,1,:),view(asm,1,:),nZ,nZ) 
    out                = Assemblystudy_scale(Lz,Lzz)
    return out,asm,dofgr
end
function zero!(out::Assemblystudy_scale)
    zero!(out.Lz )
    zero!(out.Lzz )
end
function add!(out1::Assemblystudy_scale,out2::Assemblystudy_scale) 
    add!(out1.Lz,out2.Lz)
    add!(out1.Lzz,out2.Lzz)
end
function addin!(out::Assemblystudy_scale,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxder,<:SVector{Nx}},
                                         U::NTuple{Nuder,<:SVector{Nu}},A::SVector{Na},t,SP,dbg) where{E,Nxder,Nx,Nuder,Nu,Na} # TODO make Nx,Nu,Na types
    Nz              = 2Nx+Nu+Na                        # Z =[Λ;X;U;A]       
    scaleZ          = SVector(scale.Λ...,scale.X...,scale.U...,scale.A...)
    ΔZ              = variate{2,Nz}(δ{1,Nz,𝕣}(scaleZ),scaleZ)                 
    iλ,ix,iu,ia     = gradientpartition(Nx,Nx,Nu,Na) # index into element vectors ΔZ and Lz
    ΔΛ,ΔX,ΔU,ΔA     = view(ΔZ,iλ),view(ΔZ,ix),view(ΔZ,iu),view(ΔZ,ia) 
    L,FB            = getlagrangian(eleobj, ∂0(Λ)+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA,t,SP,dbg)
    ∇L              = ∂{2,Nz}(L)
    add_value!(out.Lz ,asm[1],iele,∇L)
    add_∂!{1}( out.Lzz,asm[2],iele,∇L)
end

#------------------------------------
# an interger log10 of a number
function magnitude(x) 
    if x==0 
        NaN
    elseif abs(x)==Inf
        99
    else
        round(𝕫,log10(abs(x)))
    end
end
function short(X,n) # but 186*0.0001 = 0.018600000000000002 ...
    o   = exp10(floor(log10(X))-n+1)
    return round(Int64,X/o)*o
end
function listdoftypes(dis) # specalised for allΛXUAdofs, should be rewriten to take dofgr as input.
    type = vcat([(:Λ,f) for f∈dis.fieldX],[(:X,f) for f∈dis.fieldX],[(:U,f) for f∈dis.fieldU],[(:A,f) for f∈dis.fieldA])
    return type,unique(type)
end
# infinity norm for each dof class-type
function ∞norm(M::SparseMatrixCSC,type,types)
    ntype         = length(types)
    f             = zeros(ntype,ntype)
    for j         = 1:size(M,2)
        jtype     = findfirst(types.==type[j:j])
        for inz ∈ M.colptr[j]:M.colptr[j+1]-1
            i     = M.rowval[inz]
            itype = findfirst(types.==type[i:i])
            f[itype,jtype] = max(f[itype,jtype],abs(M.nzval[inz]))
        end
    end 
    return f
end
function ∞norm(V::Vector,type,types)    
    ntype         = length(types)
    f             = zeros(ntype)
    for i         = 1:length(V)
        itype     = findfirst(types.==type[i:i])
        f[itype]  = max(f[itype],abs(V[i]))
    end
    return f
end
"""
    scale = Muscade.study_scale(state;[SP=nothing],[verbose=false],[dbg=(;)])

Returns a named tuple of named tuples for scaling the model, accessed as
    `scaled.myclass.myfield`, for example `scale.X.tx1`.

!!! info    
    The format of `scale` is not identical to the input expected by `setscale!`

If `verbose=true`, prints out a report of the analysis underlying the proposed `scale`.  The proposed scaling depends
on the `state` passed as input - as it is computed for a given incremental matrix.
    
See also: [`setscale!`](@ref)
"""
function study_scale(state::State;SP=nothing,verbose::𝕓=true,dbg=(;))
    model,dis          = state.model,state.dis
    tmp                = state.SP 
    state.SP           = SP
    out,asm,dofgr      = prepare(Assemblystudy_scale,model,dis)
    assemble!(out,asm,dis,model,state,(dbg...,solver=:study_scale))
    state.SP           = tmp
    Z                  = zeros(getndof(dofgr))
    getdof!(state,0,Z,dofgr) 
    type,types         = listdoftypes(dis)
    matfrob            = ∞norm(out.Lzz,type,types)
    vecfrob            = ∞norm(out.Lz ,type,types)
    Zfrob              = ∞norm(Z      ,type,types)
    ntype              = length(types)
    nnz,n              = sum(matfrob.>0),length(vecfrob)
    M                  = zeros(nnz,n)
    V                  = Vector{𝕣}(undef,nnz)
    inz                = 0
    for i=1:n, j=1:n
        if matfrob[i,j]>0
            inz += 1
            V[inz]     = log10(matfrob[i,j])
            M[inz,i]   = 1
            M[inz,j]   = 1
        end
    end
    s = -(M'*M)\(M'*V)  
    S = exp10.(s)
    scaledfrob = diagm(S)*matfrob*diagm(S)
    Ss = short.(S,2)


    Xtypes = unique(dis.fieldX)
    Utypes = unique(dis.fieldU)
    Atypes = unique(dis.fieldA)
    nX     = length(Xtypes) 
    nU     = length(Utypes) 
    nA     = length(Atypes) 
    scaleΛ = (; zip(Xtypes, Ss[1:nX])...)
    scaleX = (; zip(Xtypes, Ss[nX+1:2nX])...)
    scaleU = (; zip(Utypes, Ss[2nX+1:2nX+nU])...)
    scaleA = (; zip(Atypes, Ss[2nX+nU+1:2nX+nU+nA])...)
    scale  = (Λ=scaleΛ,X=scaleX,U=scaleU,A=scaleA)

    if verbose       
        @printf "\nlog₁₀ of the ∞-norms of the blocks of the Hessian (as computed with the current scaling of the model):\n\n                    "
        for jtype = 1:ntype
            @printf "%1s%-4s " types[jtype][1] types[jtype][2]
        end
        @printf "\n"
        for itype = 1:ntype
            @printf "    %2s-%-8s  " types[itype][1] types[itype][2]
            for jtype = 1:ntype
                if matfrob[itype,jtype]==0
                    @printf "     ⋅"
                else
                    @printf "%6i" magnitude(matfrob[itype,jtype])
                end
            end
            @printf "\n"
        end
        @printf "\nlog₁₀ of the ∞-norms of the blocks of the gradient (as computed with the current scaling of the model):\n\n                    "
        for jtype = 1:ntype
            @printf "%1s%-4s " types[jtype][1] types[jtype][2]
        end
        @printf "\n                 "
        for itype = 1:ntype
            if vecfrob[itype]==0
                @printf "     ."
            else
                @printf "%6i" magnitude(vecfrob[itype])
            end
        end
        @printf "\n\nlog₁₀ of the ∞-norms of the blocks of the dofs (as computed with the current scaling of the model):\n\n                    "
        for jtype = 1:ntype
            @printf "%1s%-4s " types[jtype][1] types[jtype][2]
        end
        @printf "\n                 "
        for itype = 1:ntype
            if vecfrob[itype]==0
                @printf "     ."
            else
                @printf "%6i" magnitude(Zfrob[itype])
            end
        end
        @printf "\n\nMagnitudes of the scaling:\n\n                    "
        for jtype = 1:ntype
            @printf "%1s%-4s " types[jtype][1] types[jtype][2]
        end
        @printf "\n                 "
        for itype = 1:ntype
             @printf "%6i" s[itype]
        end
        @printf "\n\nlog₁₀ of the ∞-norms of the blocks of the SCALED Hessian:\n\n                    "
        for jtype = 1:ntype
            @printf "%1s%-4s " types[jtype][1] types[jtype][2]
        end
        @printf "\n"
        for itype = 1:ntype
            @printf "    %2s-%-8s  " types[itype][1] types[itype][2]
            for jtype = 1:ntype
                if scaledfrob[itype,jtype]==0
                    @printf "     ."
                else
                    @printf "%6i" magnitude(scaledfrob[itype,jtype])
                end
            end
            @printf "\n"
        end
        @printf "\n\nlog₁₀ of the condition number of the matrix of ∞-norms of the blocks of the SCALED Hessian = %i\n" magnitude(cond(scaledfrob))
        @printf "log₁₀ of the condition number of the matrix of ∞-norms of the blocks of the Hessian = %i\n\n" magnitude(cond(matfrob))
    end    
    return scale
end


################# study_singular

mutable struct Assemblystudy_singular{Tzz,Ticlasses,Tjclasses}  <:Assembly
    Lij   :: Tzz
    iclasses :: Ticlasses
    jclasses :: Tjclasses
end   
function prepare(::Type{Assemblystudy_singular},model,dis,iclasses=(Λ,:X,:U,:A),jclasses=iclasses) 
    idofgr             = selecteddofs(model,dis,iclasses) 
    jdofgr             = selecteddofs(model,dis,jclasses) 
    ni                 = getndof(idofgr)
    nj                 = getndof(jdofgr)
    narray,neletyp     = 3,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    _Li                = asmvec!(view(asm,1,:),idofgr,dis) 
    _Lj                = asmvec!(view(asm,2,:),jdofgr,dis) 
    Lij                = asmmat!(view(asm,3,:),view(asm,1,:),view(asm,2,:),ni,nj) 
    out                = Assemblystudy_singular(Lij,iclasses,jclasses)
    return out,asm,idofgr,jdofgr
end
function zero!(out::Assemblystudy_singular)
    zero!(out.Lij )
end
function add!(out1::Assemblystudy_singular,out2::Assemblystudy_singular) 
    add!(out1.Lz,out2.Lij)
end
function addin!(out::Assemblystudy_singular,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxder,<:SVector{Nx}},
                                         U::NTuple{Nuder,<:SVector{Nu}},A::SVector{Na},t,SP,dbg) where{E,Nxder,Nx,Nuder,Nu,Na} # TODO make Nx,Nu,Na types

    Nz              = 2Nx+Nu+Na       
    scaleZ          = SVector(scale.Λ...,scale.X...,scale.U...,scale.A...)
    ΔZ              = variate{2,Nz}(δ{1,Nz,𝕣}(scaleZ),scaleZ)                 
    iλ,ix,iu,ia     = gradientpartition(Nx,Nx,Nu,Na) # index into element vectors ΔZ and Lz
    ΔΛ,ΔX,ΔU,ΔA     = view(ΔZ,iλ),view(ΔZ,ix),view(ΔZ,iu),view(ΔZ,ia)
    L,FB            = getlagrangian(eleobj, ∂0(Λ)+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA,t,nothing,SP,dbg)
    ∇L              = ∂{2,Nz}(L)
    i               = vcat(:Λ∈out.iclasses ? iλ : 𝕫[],:X∈out.iclasses ? ix : 𝕫[],:U∈out.iclasses ? iu : 𝕫[],:A∈out.iclasses ? ia : 𝕫[])
    j               = vcat(:Λ∈out.jclasses ? iλ : 𝕫[],:X∈out.jclasses ? ix : 𝕫[],:U∈out.jclasses ? iu : 𝕫[],:A∈out.jclasses ? ia : 𝕫[])
    add_∂!{1}( out.Lij,asm[3],iele,∇L,i,j)
end
"""
    matrix = Muscade.study_singular(state;SP,[iclasses=(Λ,:X,:U,:A)],[jclasses=iclasses],[verbose::𝕓=true],[dbg=(;)])

Generates an incremental matrix for `state` (no time derivatives) corresponding to the classes required, 
and report on the null space of the matrix.

In teh present implementation, the incremental matrix is converted to full format, limiting the applicability to small models.

The function returns the incremental matrix.
"""
function study_singular(state::State;SP,iclasses=(Λ,:X,:U,:A),jclasses=iclasses,verbose::𝕓=true,dbg=(;))
    model,dis              = state.model,state.dis
    out,asm,idofgr,jdofgr  = prepare(Assemblystudy_singular,model,dis,iclasses,jclasses)
    tmp = state.SP
    state.SP = SP
    assemble!(out,asm,dis,model,state,(dbg...,solver=:study_scale))
    state.SP = tmp
    kernel =  nullspace(Matrix(out.Lij))
    nker   = size(kernel,2)
    mech   = [deepcopy(state) for iker=1:nker]
    if nker==0
        @printf("\n\nThe null-space is empty\n\n")
    end    
    for iker = 1:nker
        @printf("\n\nBase vector #%i of the null-space:\n\n",iker)
        set!(mech[iker],0,view(kernel,:,iker),jdofgr)
        describe(mech[iker])
    end
    return out.Lij
end

######### plot_matrix_sparsity
"""
    Muscade.plot_matrix_sparsity(M)

Opens a `GLMakie` figure and plots the sparsity pattern of `M::SparseMatrixCSC`.

Actual non zero-elements are plotted in green.  Remaining structuraly non-zero elements
are plotted in red.

Optional inputs:
- `size=500`        Size in pizel of the figure window.
- `title=nothing`   Title of the figure window.
- `markersize=3`    Size of dots for non-zero elements.
- `atol=1e-9`        Tolerance for actual non-zero elements.

"""
function plot_matrix_sparsity(M::SparseMatrixCSC;size=500,title=nothing,markersize=3,rtol=1e-9)
    (i,j,v)  = findnz(M)
    atol  = rtol*maximum(abs,v)
    nz = findall(abs.(v).>atol)
    if title==nothing
        title = @sprintf("nnz=%i (structural),nnz=%i (actual)",nnz(M),length(nz))
    end
    fig      = Figure(;size=(size,size))
    display(fig) # open interactive window (gets closed down by "save")
    axe      = Axis(fig[1,1];title,xticksvisible=true,yticksvisible=true,yreversed=true)
    scatter!(axe,i,j;color=:red,markersize)
    scatter!(axe,i[nz],j[nz];color=:green,markersize=2*markersize)
    #hidespines!(axe)
    return fig
end


##############  plot_block_matrix_sparsity


Base.zero(::Type{<:SparseArrays.SparseMatrixCSC})=nothing
"""
    Muscade.plot_block_matrix_sparsity(M)

Specialised tool to visualise the sparsity pattern of a matrix produced by [`DirectXUA`](@ref).
`M` is either a `Matrix` or a `SparseMatrixCSC` (the block structure), whose entries 
are themselve `SparseMatrixCSC` (the structure of each block).    

Optional inputs:
- `size=500`        Size in pizel of the figure window.
- `markersize=3`    Size of dots for non-zero elements.

"""
function plot_block_matrix_sparsity(pattern::AbstractMatrix{SparseMatrixCSC{Tv,Ti}};pixels=500,markersize=3) where{Tv,Ti}
    nbr,nbc                = size(pattern)  
    # determine the number rows in each row of blocks, store in pgr
    pgr                     = Vector{Int64}(undef,nbr+1)         # pgr[ibr]→igr pointers to the start of each block in global solution vector, where global*solution=rhs
    pgr[1]                  = 1
    for ibr                 = 1:nbr
        nlr                 = 0
        for ibc             = 1:nbc
            b = pattern[ibr,ibc]
            if ~isnothing(b)
                nlr         = size(b,1)
                break
            end
        #    ibc<nbc || @printf("row %i in pattern is empty\n",ibr)
        end
        #nlr > 0 || @printf("block-row %i has zero assigned Sparse\n",ibr)
        pgr[ibr+1]          = pgr[ibr]+nlr
    end 
    ngr                     = pgr[end]-1

    # determine the number columns in each column of blocks, store in pgc
    pgc                     = Vector{Int64}(undef,nbc+1)         # pgc[ibc]→igc pointers to the start of each block in global rhs vector
    pgc[1]                  = 1
    for ibc                 = 1:nbc
        nlc                 = 0
        for ibr             = 1:nbr
            b = pattern[ibr,ibc]
            if ~isnothing(b)
                nlc         = size(b,2)
                break
            end
            #ibr<nbr || @printf("column %i in pattern is empty\n",ibc)
        end
        #nlc > 0 || @printf("block-column %i has zero assigned Sparse\n",ibc)
        pgc[ibc+1]          = pgc[ibc]+nlc
    end 
    ngc                     = pgc[end]-1

    (I,J,B) = findnz(pattern) 

    fig      = Figure(;size=(pixels,pixels))
    display(fig) # open interactive window (gets closed down by "save")
    axe      = Axis(fig[1,1];xticksvisible=true,yticksvisible=true,yreversed=true)

    ngv = 1
    for ib ∈ eachindex(I)
        clr = iseven(I[ib]-J[ib]) ? :red : :cyan
        (i,j,v) = findnz(B[ib])
        ngv += length(i)
        li = pgc[I[ib]]-1/3
        hi = pgc[I[ib]+1]-2/3
        lj = pgr[J[ib]]-1/3
        hj = pgr[J[ib]+1]-2/3
        lines!(axe,[li,li,hi,hi,li],[lj,hj,hj,lj,lj];color=:lightgrey,linewidth=1)
        scatter!(axe,i.+(pgc[I[ib]]-1),j.+(pgr[J[ib]]-1);color=:black,markersize)
    end
    return fig
end
"""
    print_nz(S::SparseMatrixCSC)

List the structuraly non-zero entries of the sparse matrix.    
"""
function print_nz(S)
    (i,j,v) = findnz(S)
    for ℓ = 1:nnz(S)
        @printf("i=%i3, j=%i3, v=%g\n",i[ℓ],j[ℓ],v[ℓ])
    end
    @printf("_______________\n")
end

##############  Monitor

"""
    Muscade.Monitor <: AbstractElement

An element for for monitoring inputs to and outputs from
another element, during an analysis.     

Instead of adding the element to be monitored directly into the model,
add this element with the element to be monitored as argument.

Inputs and outputs are @show'n. 

# Named arguments to the constructor

- `ElementType`         The the type of element to be monitored-
- `trigger`             A function that takes `dbg` as an input and returns a boolean 
                        (`true`) to printout.
- `elementkwargs`       a `NamedTuple` containing the named arguments of the `ElementType` constructor.

"""
struct Monitor{Teleobj,Ttrigger} <: AbstractElement
    eleobj   :: Teleobj
    trigger  :: Ttrigger
end
function Monitor(nod::Vector{Node};ElementType,trigger::Function,elementkwargs)
    eleobj = ElementType(nod;elementkwargs...)
    return Monitor(eleobj,trigger)
end
doflist( ::Type{<:Monitor{Teleobj}}) where{Teleobj} = doflist(Teleobj)
@espy function lagrangian(o::Monitor{Teleobj}, Λ,X,U,A,t,SP,dbg)  where{Teleobj}
    L,FB = getlagrangian(o.eleobj,Λ,X,U,A,t,SP,(dbg...,via=Monitor))
    if o.trigger(dbg)
        @show dbg
        @show SP
        @show VALUE(Λ)
        @show VALUE(X[1])
        @show VALUE(U[1])
        @show VALUE(A)
        @show Teleobj
        @show doflist(Teleobj)
        @show L
    end
    return L,FB
end

############## @typeof

"""
    inftyp,rettyp = Muscade.@typeof(foo(args...[;kwargs...]))

    Determine the inferred type and the returned type of the output[s] returned by the relevant method-instance of foo.
    Useful to study type-stability in `lagrangian`, `residual` and more.
    This does not work on MacOS, and should thus only be used for debugging.  In tests, use `Test.@inferred`
    
"""
macro typeof(ex) # see Base.return_types, and Test.jl/@inferred
    _inferred_type(ex, __module__)
end
function _inferred_type(ex, mod)
    if Meta.isexpr(ex, :ref)
        ex = Expr(:call, :getindex, ex.args...)
    end
    Meta.isexpr(ex, :call)|| error("@inferred requires a call expression")
    farg = ex.args[1]
    if isa(farg, Symbol) && farg !== :.. && first(string(farg)) == '.'
        farg = Symbol(string(farg)[2:end])
        ex = Expr(:call, GlobalRef(Test, :_materialize_broadcasted),
            farg, ex.args[2:end]...)
    end
    result = let ex = ex
        quote
            $(if any(@nospecialize(a)->(Meta.isexpr(a, :kw) || Meta.isexpr(a, :parameters)), ex.args)
                # Has keywords
                args   = gensym()
                kwargs = gensym()
                quote
                    $(esc(args)), $(esc(kwargs)), result = $(esc(Expr(:call, _args_and_call, ex.args[2:end]..., ex.args[1])))
                    inftype = $(gen_call_with_extracted_types(mod, Base.infer_return_type, :($(ex.args[1])($(args)...; $(kwargs)...)); is_source_reflection = false))
                end
            else
                # No keywords
                quote
                    args    = ($([esc(ex.args[i]) for i = 2:length(ex.args)]...),)
                    result  = $(esc(ex.args[1]))(args...)
                    inftype = Base.infer_return_type($(esc(ex.args[1])), Base.typesof(args...))
                end
            end)
            rettype = result isa Type ? Type{result} : typeof(result)
            (inftype,rettype)
        end
    end
    return result
end

############## print_element_array

"""
    Muscade.print_element_array(eleobj,class,V)

Show a vector (or a matrix) `V`, the rows of `V` being described as corresponding to `eleobj` dof of class `class` (`:X`, `:U` or `:A`).
This can be used to print degrees of freedom, residuals, their derivatives, or gradients and Hessian of the Lagrangian.

See also: [`diffed_residual`](@ref), [`diffed_lagrangian`](@ref)
"""    
print_element_array(ele::AbstractElement,class::Symbol,V::AbstractVector) = print_element_array(ele,class,reshape(V,(length(V),1)))
function print_element_array(ele::Eletyp,class::Symbol,V::AbstractMatrix) where{Eletyp<:AbstractElement}
    inod,~,field     = Muscade.getdoflist(Eletyp)
    iVdof            = Muscade.getidof(Eletyp,class)
    (nV,ncol)      = size(V)
    @assert nV==Muscade.getndof(Eletyp,class)
    @printf "    i  ieldof               doftyp   inod |"
    for icol = 1:ncol
        @printf "  %10i" icol 
    end
    @printf "\n__________________________________________|"
    for icol = 1:ncol
        @printf "____________" 
    end
    @printf "\n"

    for iV ∈ 1:nV
        idof = iVdof[iV]
        @printf " %4d    %4d     %16s  %5d |" iV idof field[idof] inod[idof] 
        for icol = 1:ncol
            @printf "  %10.3g" V[iV,icol] 
        end
        @printf "\n"
    end
end

############ diffed_lagrangian

"""
    Muscade.diffed_lagrangian{P}(eleobj;Λ,X,U,A,t=0.,SP=nothing)

Compute the Lagrangian, its gradients and Hessian, and the memory of an element.
For element debugging and testing. 

`P`, the order of differentiation must be 1 or 2.

The output is a `NamedTuple` with fields `Λ`, `X`, `U`, `A`, `t`, `SP` echoing the inputs and fields
- `∇L` of format `∇L[iclass][ider]`so that for example `∇L[2][3]` contains the gradient of the Lagrangian wrt to the acceleration.
   `iclass` is 1,2,3 and 4 for `Λ`, `X`, `U` and `A` respectively.
- if `P==2`: `HL` of format `HL[iclass,jclass][ider,jder]`so that for example `HL[1,2][1,3]` contains the mass matrix.
- `FB` as returned by `lagrangian`

See also: [`diffed_residual`](@ref), [`print_element_array`](@ref)
"""     
struct diffed_lagrangian{P} end
function diffed_lagrangian{P}(ele::Eletyp; Λ,X,U,A, t::𝕣=0.,SP=nothing) where{P,Eletyp<:AbstractElement}
    OX,OU,IA         = length(X)-1,length(U)-1,1
    Nλ               = length(   Λ ) 
    Nx               = length(∂0(X)) 
    Nu               = length(∂0(U)) 
    Na               = length(   A ) 

    if (Nλ,Nx,Nu,Na) ≠ getndof(Eletyp,(:X,:X,:U,:A))
        display(Eletyp)
        @printf("diffed_lagrangian received %i Λ, %i X, %i U and %i A dofs\n",Nλ,Nx,Nu,Na)
        @printf("element requires           %i Λ, %i X, %i U and %i A dofs\n",getndof(Eletyp,:X),getndof(Eletyp,:X),getndof(Eletyp,:U),getndof(Eletyp,:A))
        @assert false
    end
    λxua      = ( 1,    2,    3,  4)
    ndof      = (Nx,   Nx,   Nu, Na)
    nder      = ( 1, OX+1, OU+1, IA)
    Np        = Nx + Nx*(OX+1) + Nu*(OU+1) + Na*IA # number of partials 
    d         = revariate{P}((Λ=Λ,X=X,U=U,A=A))

    L,FB      = lagrangian(ele, d.Λ,d.X,d.U,d.A,t,SP,(;calledby=:test_element))
    
    if P==1
        ∇Lz       = ∂{1,Np}(L)
        ∇L        = Vector{Vector{Any}}(undef,4  )
        pα        = 0   # points into the partials, 1 entry before the start of relevant partial derivative in α,ider-loop
        for α∈λxua 
            ∇L[α] = Vector{Any}(undef,nder[α])
            for i=1:nder[α]   # we must loop over all time derivatives to correctly point into the adiff-partials...
                iα       = pα.+(1:ndof[α])
                pα      += ndof[α]
                ∇L[α][i] = ∇Lz[iα]
                pβ       = 0
            end
        end
        return (Λ=Λ,X=X,U=U,A=A,t=t,SP=SP,∇L=∇L,FB=FB)#,inftyp=inftyp,rettyp=rettyp)
    elseif P==2
        ∇Lz,HLz   = value_∂{1,Np}(∂{2,Np}(L))
        ∇L        = Vector{Vector{Any}}(undef,4  )
        HL        = Matrix{Matrix{Any}}(undef,4,4)
        pα        = 0   # points into the partials, 1 entry before the start of relevant partial derivative in α,ider-loop
        for α∈λxua 
            ∇L[α] = Vector{Any}(undef,nder[α])
            for i=1:nder[α]   # we must loop over all time derivatives to correctly point into the adiff-partials...
                iα       = pα.+(1:ndof[α])
                pα      += ndof[α]
                ∇L[α][i] = ∇Lz[iα]
                pβ       = 0
                for β∈λxua 
                    HL[α,β] = Matrix{Any}(undef,nder[α],nder[β])
                    for j=1:nder[β]
                        iβ   = pβ.+(1:ndof[β])
                        pβ  += ndof[β]
                        HL[α,β][i,j] = HLz[iα,iβ]
                    end
                end
            end
        end
        return (Λ=Λ,X=X,U=U,A=A,t=t,SP=SP,∇L=∇L,HL=HL,FB=FB)#,inftyp=inftyp,rettyp=rettyp)
    end
end

############# diffed_residual 

"""
    Muscade.diffed_residual(eleobj;X,U,A,t=0.,SP=nothing)

Compute the residual, its gradients, and the memory of an element.
For element debugging and testing. 

The output is a `NamedTuple` with fields `X`, `U`, `A`, `t`, `SP` echoing the inputs and fields
- `R` containing the residual
- `∇R` of format `∇R[iclass][ider]`so that for example `∇R[2][3]` contains the mass matrix.
   `iclass` is 2,3 and 4 for `X`, `U` and `A` respectively.
- `FB` as returned by `residual`

See also: [`diffed_lagrangian`](@ref), [`print_element_array`](@ref)
"""     
function diffed_residual(ele::Eletyp; X,U,A, t::𝕣=0.,SP=nothing) where{Eletyp<:AbstractElement}
    OX,OU,IA         = length(X)-1,length(U)-1,1
    Nx               = length(∂0(X)) 
    Nu               = length(∂0(U)) 
    Na               = length(   A ) 
    if (Nx,Nu,Na) ≠ getndof(Eletyp,(:X,:U,:A))
        display(Eletyp)
        @printf("diffed_residual received %i X, %i U and %i A dofs\n",Nx,Nu,Na)
        @printf("element requires         %i X, %i U and %i A dofs\n",getndof(Eletyp,:X),getndof(Eletyp,:U),getndof(Eletyp,:A))
        @assert false
    end

    xua       = (    2,    3,  4)
    ndof      = (0, Nx,   Nu, Na)
    nder      = (0 ,OX+1, OU+1, IA)
    Np        = Nx*(OX+1) + Nu*(OU+1) + Na*IA # number of partials 
    d         = revariate{1}((X=X,U=U,A=A))
    r_,FB     = residual(ele, d.X,d.U,d.A,t,SP,(;calledby=:test_element))
    #inftyp,rettyp = @typeof(residual(ele, X∂,U∂,A∂,t,SP,(;calledby=:test_element)))
    R,∇r      = value_∂{1,Np}(r_)

    ∇R        = Vector{Vector{Any}}(undef,4  )
    pα        = 0   # points into the partials, 1 entry before the start of relevant partial derivative in α,ider-loop
    for α∈xua 
        ∇R[α] = Vector{Any}(undef,nder[α])
        for i=1:nder[α]   # we must loop over all time derivatives to correctly point into the adiff-partials...
            iα       = pα.+(1:ndof[α])
            pα      += ndof[α]
            ∇R[α][i] = ∇r[:,iα]
        end
    end
    return (X=X,U=U,A=A,t=t,SP=SP,R=R,∇R=∇R,FB=FB)#,inftyp=inftyp,rettyp=rettyp)
end

