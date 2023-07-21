
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
    ΔΛ,ΔX,ΔU,ΔA     = view(ΔZ,iλ),view(ΔZ,ix),view(ΔZ,iu),view(ΔZ,ia) # TODO Static?
    L,χn,FB         = getlagrangian(implemented(eleobj)...,eleobj, Λ+ΔΛ, (∂0(X)+ΔX,),(∂0(U)+ΔU,),A+ΔA,t,nothing,nothing,SP,dbg)
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
function frobenius(M::SparseMatrixCSC,type,types)
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
function frobenius(V::Vector,type,types)    
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
function studyscale(state::State;verbose::𝕓=true,dbg=(;))
    model,dis          = state.model,state.dis
    out,asm,dofgr      = prepare(AssemblyStudyScale,model,dis)
    assemble!(out,asm,dis,model,state,(dbg...,solver=:studyscale))

    type,types         = listdoftypes(dis)
    matfrob            = frobenius(out.Lzz,type,types)
    vecfrob            = frobenius(out.Lz ,type,types)
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
        @printf "\nMagnitudes of the Frobenius norms of the blocks of the Hessian (as computed with the current scaling of the model):\n\n                   "
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
        @printf "\nMagnitudes of the Frobenius norms of the blocks of the gradient(as computed with the current scaling of the model):\n\n                 "
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
        @printf "\n\nMagnitude of the condition number of the matrix of Frobenius norms of the blocks of the Hessian = %i\n\n" magnitude(cond(matfrob))
        @printf "\nMagnitudes of the Frobenius norms of the blocks of the SCALED Hessian:\n\n                 "
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
        @printf "\n\nMagnitude of the condition number of the matrix of Frobenius norms of the blocks of the SCALED Hessian = %i\n\n" magnitude(cond(scaledfrob))
    end    
    return scale
end


