module TestEspy
using Muscade,Test#,StaticArrays

## Test request

struct Element  end
requestable(el::Element)= (X=(3,2),  gp=forloop(2, (F=(2,2), material=(ε=(2,2),Σ=(2,2)))) )
el = Element()
r8 = @request X,gp[].(F,material.(Σ,ε))

ra    = requestable(el)

@testset "requestable" begin
    @test ra.X == (3,2)
    @test ra.gp.body == (F = (2, 2), material = (ε = (2, 2), Σ = (2, 2)))
end

@testset "Makekey" begin
    @test Muscade.makekey_symbol(2,(2,3)) == ([3 5 7; 4 6 8], 8)
    @test Muscade.makekey_tuple(2,:((a,b)),(a=(2,3),b=(2,2))) == ((a=[3 5 7; 4 6 8], b=[9 11; 10 12]), 12)
    @test Muscade.makekey_tuple(2,:((a,b)),(a=(2,3),b=scalar)) == ((a=[3 5 7; 4 6 8], b=9), 9)
    @test Muscade.makekey_loop(2,:((a,b)),forloop(3,(a=(2,3),b=(2,2)))) == ([(a=[3 5 7; 4 6 8], b=[9 11; 10 12]), (a=[13 15 17; 14 16 18], b=[19 21; 20 22]), (a=[23 25 27; 24 26 28], b=[29 31; 30 32])], 32)
    @test Muscade.makekey_tuple(2,:((a[].(b,)),),(a=forloop(2,(b=(2,3),)),)) == (( a=[(b=[3 5 7; 4 6 8],), (b=[9 11 13; 10 12 14],)] ,), 14)
    @test Muscade.makekey_tuple(2,:(a[].b),(a=forloop(2,(b=(2,3),)),)) == ((a=[(b=[3 5 7; 4 6 8],), (b=[9 11 13; 10 12 14],)],), 14)
    @test makekey(r8,ra) == ((X = [1 4; 2 5; 3 6], gp = [(F = [7 9; 8 10], material = (Σ = [11 13; 12 14], ε = [15 17; 16 18])), (F = [19 21; 20 22], material = (Σ = [23 25; 24 26], ε = [27 29; 28 30]))]), 30)
end

## annotated code (written by user)

struct El  end
requestable(el::El)= (gp=forloop(2, (z=scalar,s=scalar, material=(a=scalar,b=scalar))),)
el = El()


@espy function residual(x,y)
    ngp=2
    r = 0
    for igp=1:ngp
        :z = x[igp]+y[igp]
        :s,dum  = :material(z)
        r += s
    end
    return r
end
@espy function material(z)
    :a = z+1
    :b = a*z
    return b,3.
end

## Test result extraction
req       = @request gp[].(s,z,material.(a,b))  # generates an expression, not a variable
key,nkey  = makekey(req,requestable(el))
nstep,nel = 2,3
out       = Array{Float64,3}(undef,nkey,nel,nstep)
istep,iel = 1,2
x,y       = [1.,2.],[.5,.2]
r         = residual(@view(out[:,iel,istep]),key,x,y)

@testset "Espy" begin
    @test key == (gp=[(s=1, z=2, material=(a=3, b=4)), (s=5, z=6, material=(a=7, b=8))],)
    @test out[:,iel,istep] ≈ [3.75, 1.5, 2.5, 3.75, 7.04, 2.2, 3.2, 7.04]
end
out = Vector{Float64}(undef,5)
out[[1,2]] .= [1,2]
out[[3,4]] .= (3,4)
# out[[5,6]] .= @SVector [5,6] works, but requires StaticArrays in the test environment
out[[5]]   .= 5
@testset "Assign to out" begin
    @test out == [1.,2.,3.,4.,5.]
end

end
