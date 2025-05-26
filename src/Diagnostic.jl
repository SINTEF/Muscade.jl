
### getting printout of the model
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

####################

mutable struct AssemblyStudyScale{Tz,Tzz}  <:Assembly
    Lz    :: Tz
    Lzz   :: Tzz
end   
function prepare(::Type{AssemblyStudyScale},model,dis) 
    dofgr              = allΛXUAdofs(model,dis)
    nZ                 = getndof(dofgr)
    narray,neletyp     = 2,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    Lz                 = asmvec!(view(asm,1,:),dofgr,dis) 
    Lzz                = asmmat!(view(asm,2,:),view(asm,1,:),view(asm,1,:),nZ,nZ) 
    out                = AssemblyStudyScale(Lz,Lzz)
    return out,asm,dofgr
end
function zero!(out::AssemblyStudyScale)
    zero!(out.Lz )
    zero!(out.Lzz )
end
function add!(out1::AssemblyStudyScale,out2::AssemblyStudyScale) 
    add!(out1.Lz,out2.Lz)
    add!(out1.Lzz,out2.Lzz)
end
function addin!(out::AssemblyStudyScale,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxder,<:SVector{Nx}},
                                         U::NTuple{Nuder,<:SVector{Nu}},A::SVector{Na},t,SP,dbg) where{E,Nxder,Nx,Nuder,Nu,Na} # TODO make Nx,Nu,Na types
    Nz              = 2Nx+Nu+Na                        # Z =[Λ;X;U;A]       
    scaleZ          = SVector(scale.Λ...,scale.X...,scale.U...,scale.A...)
    ΔZ              = variate{2,Nz}(δ{1,Nz,𝕣}(scaleZ),scaleZ)                 
    iλ,ix,iu,ia     = gradientpartition(Nx,Nx,Nu,Na) # index into element vectors ΔZ and Lz
    ΔΛ,ΔX,ΔU,ΔA     = view(ΔZ,iλ),view(ΔZ,ix),view(ΔZ,iu),view(ΔZ,ia) 
    L,FB         = getlagrangian(eleobj, ∂0(Λ)+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA,t,SP,dbg)
    ∇L              = ∂{2,Nz}(L)
    add_value!(out.Lz ,asm[1],iele,∇L)
    add_∂!{1}( out.Lzz,asm[2],iele,∇L)
end

#------------------------------------
magnitude(x) = x==0 ? NaN : round(𝕫,log10(abs(x)))
function short(X,n) # but 186*0.0001 = 0.018600000000000002 ...
    o   = exp10(floor(log10(X))-n+1)
    return round(Int64,X/o)*o
end
function listdoftypes(dis) # specalised for allΛXUAdofs, should be rewriten to take dofgr as input.
    type = vcat([(:Λ,f) for f∈dis.fieldX],[(:X,f) for f∈dis.fieldX],[(:U,f) for f∈dis.fieldU],[(:A,f) for f∈dis.fieldA])
    return type,unique(type)
end
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
    scale = studyscale(state;[verbose=false],[dbg=(;)])

Returns a named tuple of named tuples for scaling the model, accessed as
    `scaled.myclass.myfield`, for example `scale.X.tx1`.

!!! info    
    Currently, the format of `scale` is not identical to the input expected by `setscale!`: work in progress

If `verbose=true`, prints out a report of the analysis underlying the proposed `scale`.  The proposed scaling depends
on the `state` passed as input - as it is computed for a given incremental matrix.
    
See also: [`setscale!`](@ref)
"""
function studyscale(state::State;SP,verbose::𝕓=true,dbg=(;))
    model,dis          = state.model,state.dis
    tmp                = state.SP 
    state.SP           = SP
    out,asm,dofgr      = prepare(AssemblyStudyScale,model,dis)
    assemble!(out,asm,dis,model,state,(dbg...,solver=:studyscale))
    state.SP           = tmp
    type,types         = listdoftypes(dis)
    matfrob            = ∞norm(out.Lzz,type,types)
    vecfrob            = ∞norm(out.Lz ,type,types)
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
        @printf "\nlog₁₀ of the ∞-norms of the blocks of the Hessian (as computed with the current scaling of the model):\n\n                   "
        for jtype = 1:ntype
            @printf "%1s%-4s " types[jtype][1] types[jtype][2]
        end
        @printf "\n"
        for itype = 1:ntype
            @printf "    %2s-%-8s  " types[itype][1] types[itype][2]
            for jtype = 1:ntype
                if matfrob[itype,jtype]==0
                    @printf "     ."
                else
                    @printf "%6i" magnitude(matfrob[itype,jtype])
                end
            end
            @printf "\n"
        end
        @printf "\nlog₁₀ of the ∞-norms of the blocks of the gradient(as computed with the current scaling of the model):\n\n                 "
        for itype = 1:ntype
            if vecfrob[itype]==0
                @printf "     ."
            else
                @printf "%6i" magnitude(vecfrob[itype])
            end
        end
        @printf "\n\nMagnitudes of the scaling:\n\n                 "
        for itype = 1:ntype
             @printf "%6i" s[itype]
        end
        @printf "\n\nlog₁₀ of the condition number of the matrix of ∞-norms of the blocks of the Hessian = %i\n\n" magnitude(cond(matfrob))
        @printf "\nlog₁₀ of the ∞-norms of the blocks of the SCALED Hessian:\n\n                 "
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
        @printf "\n\nlog₁₀ of the condition number of the matrix of ∞-norms of the blocks of the SCALED Hessian = %i\n\n" magnitude(cond(scaledfrob))
    end    
    return scale
end


#################

mutable struct AssemblyStudySingular{Tzz,Ticlasses,Tjclasses}  <:Assembly
    Lij   :: Tzz
    iclasses :: Ticlasses
    jclasses :: Tjclasses
end   
function prepare(::Type{AssemblyStudySingular},model,dis,iclasses=(Λ,:X,:U,:A),jclasses=iclasses) 
    idofgr             = selecteddofs(model,dis,iclasses) 
    jdofgr             = selecteddofs(model,dis,jclasses) 
    ni                 = getndof(idofgr)
    nj                 = getndof(jdofgr)
    narray,neletyp     = 3,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    _Li                = asmvec!(view(asm,1,:),idofgr,dis) 
    _Lj                = asmvec!(view(asm,2,:),jdofgr,dis) 
    Lij                = asmmat!(view(asm,3,:),view(asm,1,:),view(asm,2,:),ni,nj) 
    out                = AssemblyStudySingular(Lij,iclasses,jclasses)
    return out,asm,idofgr,jdofgr
end
function zero!(out::AssemblyStudySingular)
    zero!(out.Lij )
end
function add!(out1::AssemblyStudySingular,out2::AssemblyStudySingular) 
    add!(out1.Lz,out2.Lij)
end
function addin!(out::AssemblyStudySingular,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxder,<:SVector{Nx}},
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
    matrix = studysingular(state;SP,[iclasses=(Λ,:X,:U,:A)],[jclasses=iclasses],[verbose::𝕓=true],[dbg=(;)])

Generates an incremental matrix for `state` (no time derivatives) corresponding to the classes required, 
and report on the null space of the matrix.

To do so, the incremental matrix is converted to full format, limiting the applicability to small models.

The function returns the incremental matrix.
"""
function studysingular(state::State;SP,iclasses=(Λ,:X,:U,:A),jclasses=iclasses,verbose::𝕓=true,dbg=(;))
    model,dis              = state.model,state.dis
    out,asm,idofgr,jdofgr  = prepare(AssemblyStudySingular,model,dis,iclasses,jclasses)
    tmp = state.SP
    state.SP = SP
    assemble!(out,asm,dis,model,state,(dbg...,solver=:studyscale))
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

############## describe state to the user
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
"""
    describe(state;class=:all)

Provide a description of the dofs stored in `state`.
`class` can be either `:all`, `:Λ`, `:ΛX`, `:X`, `:U`, `:A` or `:scale`

See also: [`solve`](@ref)
"""
function describe(state::State;class::Symbol=:all)
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

