using Muscade, StaticArrays

#Create two functions
f(ξ)    =   SVector(ξ[1]^ξ[2],ξ[1]^2)
g(ξ)    =   SVector(log(ξ[1]),log(ξ[2]))	

# Point at which we evaluate the functions and gradients
x       =   SVector(3.,2.)
@show typeof(x)
# returns typeof(x) = SVector{2, Float64}
 
# Create an adiff xₐ at that point x
N       =   length(x)
xₐ      =   Muscade.variate_{1,N}(x)
@show typeof(xₐ)
# returns typeof(xₐ) = SVector{2, ∂ℝ{1, 2, Float64}}

# Compute the gradient of f at x 
fxₐ     =   f(xₐ)
@show value{1}(fxₐ)
# returns value{1}(fxₐ) = [9.0, 9.0]
@show ∂{1,N}(fxₐ)   
# returns ∂{1, N}(fxₐ) = [6.0 9.887510598012987; 6.0 0.0]

# Cmpute the gradient of g at x
gxₐ     =   g(xₐ)
@show value{1}(gxₐ)
# returns value{1}(gxₐ) = [1.0986122886681098, 0.6931471805599453]
@show ∂{1,N}(gxₐ)   
# returns ∂{1, N}(gxₐ) = [0.3333333333333333 0.0; 0.0 0.5]

# Cmpute the gradient of h=g∘f at x
hxₐ     =   g(f(xₐ))
@show value{1}(hxₐ)
# returns  value{1}(hxₐ) = [2.1972245773362196, 2.1972245773362196]
@show ∂{1,N}(hxₐ)   
# returns  ∂{1, N}(hxₐ) = [0.6666666666666666 1.0986122886681098; 0.6666666666666666 0.0]


# Compute the gradient of g, for a variable y depends on x
# Create an adiff yₐ at point xₐ
yₐ      =   Muscade.variate_{2,N}(fxₐ)
# Note: The level of preceedence (here, 2) can be obtained using constants(fxₐ)
@show typeof(yₐ)
# returns typeof(yₐ) = SVector{2, ∂ℝ{2, 2, ∂ℝ{1, 2, Float64}}}

gyₐ     =   g(yₐ)

# The result is an adiff contains nested partial derivatives (g wrt y and y wrt x) 
@show gyₐ_value   = value{2}(gyₐ)
# returns gyₐ_value = value{2}(gyₐ) = ∂ℝ{1, 2, Float64}[2.2+∂₁⟨0.667,1.1⟩ , 2.2+∂₁⟨0.667,0⟩ ]
# The gradient of g at yₐ also contains information about the gradient of y wrt x
@show ∂{2,N}(gyₐ);
# returns ∂{2, N}(gyₐ) = ∂ℝ{1, 2, Float64}[0.111+∂₁⟨-0.0741,-0.122⟩  0+∂₁⟨0,0⟩ ; 0+∂₁⟨0,0⟩  0.111+∂₁⟨-0.0741,0⟩ ]
 
# Continue to unfold (decreasing the preceedence to 1) to get the the derivative of g wrt x
@show value{1}(gyₐ_value)
# returns value{1}(gyₐ_value) = [2.1972245773362196, 2.1972245773362196]
@show ∂{1,N}(gyₐ_value);
# returns ∂{1, N}(gyₐ_value) = [0.6666666666666666 1.0986122886681098; 0.6666666666666666 0.0]
