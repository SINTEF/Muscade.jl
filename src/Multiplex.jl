using SpecialFunctions
struct Multiplex{A} <:ℝ where{A<:AbstractVector{R}} where{R<:ℝ}  
    x  :: A
end

demultiplex(a::Multiplex)     = a.x
demultiplex(a::SVector{L1,   Multiplex{SVector{M,R}}}) where{L1   ,M,R} = SMatrix{L1,M,R,L1*M}(aᵢ.x[j] for aᵢ∈a, j=1:length(a[1].x))
#demultiplex(a::SArray{S,Multiplex{SVector{M,R}}}) where{S<:Tuple{s},M,R} where{s} = SArray{Tuple{s...,M},R}(aᵢ.x[j] for aᵢ∈a, j=1:length(a[1].x))
demultiplex(a::AbstractArray) = [ aᵢ.x[j] for aᵢ∈a, j=1:length(a[1].x)]

for OP =(:zero,:one,:typemax,:typemin,:floatmax,:floatmin,:eps)
    eval(:(@inline Base.$OP(::Multiplex{A}) where{A<:AbstractVector{R}} where{R<:ℝ}  = $OP(R)))
end
for OP = (:atan,:hypot,:(+),:(-),:(*),:(/),:(^))
    eval(quote
        @inline Base.$OP(a::Multiplex,b::Multiplex) = Multiplex($OP.(a.x,b.x))
        @inline Base.$OP(a::Multiplex,b::ℝ        ) = Multiplex($OP.(a.x,b  ))
        @inline Base.$OP(a::ℝ        ,b::Multiplex) = Multiplex($OP.(a  ,b.x))
    end)
end
for OP =(:(+),:(-),:abs,:conj,:sqrt,:cbrt,:abs2,:inv,:log,:log10,:log2,:log1p,:exp,:exp2,:exp10,:expm1,:sin,:cos,:tan,:sinpi,:cospi,:sec,:csc,:cot,
         :sind,:cosd,:tand,:secd,:acsc,:acot,:asind,:acosd,:atand,:asecd,:acscd,:acotd,:sinh,:cosh,:tanh,:sech,:csch,:coth,:asinh,:acosh,:atanh,:asech,:acsch,:acoth)
    eval(:(@inline Base.$OP(a::Multiplex) = Multiplex($OP.(a.x))))
end
for OP =(:erf,:erfc,:erfi,:gamma,:lgamma,:airy,:airyprime,:airyai,:airybi,:airyaiprime,:airybiprime,:besselj0,:besselj1,:bessely0,:bessely1)
    eval(:(@inline SpecialFunctions.$OP(a::Multiplex) = Multiplex($OP.(a.x))))
end
