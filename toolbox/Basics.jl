"""
    vec3(v, ind)

Create a 3-element static vector from selected elements of vector `v`.

# Arguments
- `v`: The input vector from which elements are selected.
- `ind`: A collection of 3 indices specifying which elements of `v` to include.

# Returns
- An `SVector{3}` containing the elements `v[ind[1]]`, `v[ind[2]]`, `v[ind[3]]`.
"""
vec3(v,ind) = SVector{3}(v[i] for i∈ind);

"""
	clutch(x,x₁,x₂,y₁,y₂,γ)
	
Interpolation between x₁ and x₂ in ℝ, to (possibly vector-valued) y₁ and y₂. 
The parameter γ enables controlling how fast clutch goes from y₁ and y₂. 
A typical application of clutch are static structural anlayses, where different load components should be applied with a specific sequence/magnitude to ensure convergence of the analysis

# Inputs
- x ∈ ℝ is the point at which clutch should be evaluated
- x₁ ∈ ℝ is the point at which clutch starts. The output is equal to y₁ for x<x₁.
- x₂ ∈ ℝ is the point at which clutch stops. The output is equal to y₂ for x>x₂.
- y₁ and y₂ are the values of clutch for x<x₁ and x>x₂, respectively. 
- γ>0 controls the shape of the interpolation. 
        If γ=1, linear interpolation between x₁ and x₂ 
        If γ<1, faster progression towards y₂ first, and then slower
        If γ>1, slower progression towards y₂ first, and then faster
 
# Output
- The interpolated value at x, possibly vector-valued. 

# Example
```
x = -1:0.01:2
γ_ = [0.1,0.5,1.,2,10.]
using GLMakie
fig = Figure(size = (1000,1000))
ax = Axis(fig[1, 1],ylabel="clutch() value")
[lines!(ax,x,clutch(x,0,1,0,1,γ),label="γ="*string(γ)) for γ∈γ_]; 
axislegend(); 
xlims!(ax,-1,2); display(fig)
```
"""
function clutch(x,x₁,x₂,y₁,y₂,γ)
x₀  = min.(max.(x,x₁),x₂)
y   = y₁.+(y₂.-y₁).*((x₀.-x₁)./(x₂-x₁)).^γ
return  y 
end
