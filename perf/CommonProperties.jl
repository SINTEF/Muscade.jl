# Read RIFLEX motions, tensions, bending moment and curvatures time series (verification data)
file_path = normpath(joinpath("..","examples","SCR.csv")) 
df = CSV.read(file_path, DataFrame; delim=' ')
dynMotionX = vcat(0,df[:,"Displacement in x - direction"].-df[1,"Displacement in x - direction"],df[end,"Displacement in x - direction"]-df[1,"Displacement in x - direction"])
dynMotionZ = vcat(0,df[:,"Displacement in z - direction"].-df[1,"Displacement in z - direction"],df[end,"Displacement in z - direction"]-df[1,"Displacement in z - direction"])
dynMotionT = vcat(-10.,df[:,"x"],df[end,"x"]+100)
xMotion = linear_interpolation(dynMotionT,-dynMotionX)
zMotion = linear_interpolation(dynMotionT,dynMotionZ);
@functor with(xMotion) horizMove(x,t)= SVector( x[1]-(1.0 - clutch(t,-10.,0,0.,1.,3)*181.0 + xMotion(t)) )
@functor with(zMotion) vertMove(x,t) = SVector( x[1]-(      clutch(t,-10.,0,0.,1.,3)*303.1 + zMotion(t)) )

# Some physical constants 
const g=9.81
const ρ=1025.;

# Normalized gravity field, ramping up during static analysis
@functor with() ramping_gravityField(t)= clutch(t,-9.,-5.,0.,1.,1) * SVector(0.,0.,-1.);

# Define parameters for cross-section 1
x1_D    = 0.429                     # Outer diameter [m]
x1_t    = 0.022                     # Thickness [m]
x1_steelArea = (x1_D-x1_t)*π*x1_t   # Pipe wall area [m^2]
x1_innerArea = (x1_D-2*x1_t)^2*π/4  # Inner area [m^2]
x1_ρₛ   = 7850.                     # Steel denstiy [kg/m^3]
x1_ρᵢ   = 200.                      # Internal fluid density [kg/m^3]
x1_EA   = 5.823e9                   # Axial stiffness [N]
x1_EI₂  = 1.209e8                   # Bending stiffness [Nm/(1/m)]
x1_EI₃  = 1.209e8                   # Bending stiffness [Nm/(1/m)]
x1_GJ   = 9.347e7                   # Torsional stiffness [Nm/(rad/m)]
x1_μ    = x1_steelArea*x1_ρₛ + x1_innerArea*x1_ρᵢ   # Mass per unit length [m]
x1_ι₁   = 0.2053^2*x1_steelArea*x1_ρₛ # Moment of inertia about x-axis per unit length [kgm²/m]
x1_w    = x1_μ*g - π*x1_D^2/4*ρ*g   # Weight per unit length [N/m]
x1_Dh    = 0.459                    # Hydrodynamic outer diameter [m]
x1_Ca₂  = 1.0 *         ρ * π*x1_Dh^2/4 # Transverse added mass coefficients [N/m/(m/s^2)]
x1_Ca₃  = 1.0 *         ρ * π*x1_Dh^2/4
x1_Cq₂  = 1.0 *   0.5 * ρ * x1_Dh   # Transverse drag coefficients [N/m/(m/s)^2]
x1_Cq₃  = 1.0 *   0.5 * ρ * x1_Dh


# Define parameters for cross-section 2
x2_D    = 0.441                     # Outer diameter [m]
x2_t    = 0.028                     # Thickness [m]
x2_steelArea = (x2_D-x2_t)*π*x2_t   # Pipe wall area [m^2]
x2_innerArea = (x2_D-2*x2_t)^2*π/4  # Inner area [m^2]
x2_ρₛ   = 7850.                     # Steel denstiy [kg/m^3]
x2_ρᵢ   = 200.                      # Internal fluid density [kg/m^3]
x2_EA   = 7.520e9                   # Axial stiffness [N]
x2_EI₂  = 1.611e8                   # Bending stiffness [Nm/(1/m)]
x2_EI₃  = 1.611e8                   # Bending stiffness [Nm/(1/m)]
x2_GJ   = 1.245e8                   # Torsional stiffness [Nm/(rad/m)]
x2_μ    = x2_steelArea*x2_ρₛ + x2_innerArea*x2_ρᵢ    # Mass per unit length [m]
x2_ι₁   = 0.2084^2*x2_steelArea*x2_ρₛ # Moment of inertia about x-axis per unit length [kgm²/m]
x2_w    = x2_μ*g - π*x2_D^2/4*ρ*g   # Weight per unit length [N/m]
x2_Dh    = 0.471                    # Hydrodynamic outer diameter [m]
x2_Ca₂  = 1.0 *         ρ * π*x2_Dh^2/4 # Transverse added mass coefficients [N/m/(m/s^2)]
x2_Ca₃  = 1.0 *         ρ * π*x2_Dh^2/4
x2_Cq₂  = 1.0 *   0.5 * ρ * x2_Dh   # Transverse drag coefficients [N/m/(m/s)^2]
x2_Cq₃  = 1.0 *   0.5 * ρ * x2_Dh

# Model SCR, starting from the extremity located on seabed
nel         = [120,      60,     20] # Number of elements per segment
segLength   = [300.,    300.,   80.] # Segment lengths
