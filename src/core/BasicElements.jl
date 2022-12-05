struct XdofCost{Tcost,Field,Derivative} <: AbstractElement
    cost :: Tcost # Function 
end
XdofCost(nod::Vector{Node};field::Symbol,cost::Tcost,derivative=0::𝕫) where{Tcost<:Function} = XdofCost{Tcost,field,derivative}(cost)
Muscade.doflist(::Type{XdofCost{Tcost,Field,Derivative}}) where{Tcost,Field,Derivative} = (inod =(1,), class=(:X,), field=(Field,))
@espy function Muscade.lagrangian(o::XdofCost{Tcost,Field,Derivative}, δX,X,U,A, t,ε,dbg) where{Tcost,Field,Derivative}
        return o.cost(∂n(X,Derivative)[1])
end

struct UdofCost{Tcost,Field,Derivative} <: AbstractElement
    cost :: Tcost # Function 
end
UdofCost(nod::Vector{Node};field::Symbol,cost::Tcost,derivative=0::𝕫) where{Tcost<:Function} = UdofCost{Tcost,field,derivative}(cost)
Muscade.doflist(::Type{UdofCost{Tcost,Field,Derivative}}) where{Tcost,Field,Derivative} = (inod =(1,), class=(:U,), field=(Field,))
@espy function Muscade.lagrangian(o::UdofCost{Tcost,Field,Derivative}, δX,X,U,A, t,ε,dbg) where{Tcost,Field,Derivative}
        return o.cost(∂n(XU,Derivative)[1])
end

struct AdofCost{Tcost,Field} <: AbstractElement
    cost :: Tcost # Function 
end
AdofCost(nod::Vector{Node};field::Symbol,cost::Tcost) where{Tcost<:Function} = AdofCost{Tcost,field}(cost)
Muscade.doflist(::Type{AdofCost{Tcost,Field}}) where{Tcost,Field} = (inod =(1,), class=(:A,), field=(Field,))
@espy function Muscade.lagrangian(o::AdofCost{Tcost,Field}, δX,X,U,A, t,ε,dbg) where{Tcost,Field}
        return o.cost(A[1])
end
