"""

Muscade provides `unit` to transform quantities to and from basic SI units.
This allows to develop Muscade-based applications with a coherent unit system.

```
using Muscade
using Muscade: m, kg, pound, foot
rho          = 3←pound/foot^3                        # convert to SI
vieuxquintal = 1000*pound                            # define new unit
printf("Density [pound/foot^3] %f",rho→pound/foot^3) # convert from SI
```
Whole arrays can be converted in the same way.

In the above, `rho` would be of type `Float64`.  This contrasts with
`Unitful.jl` which would make ``rho` be of a type containing data about 
dimensionality thus providing check against operations like `length+surface`.
The drawback is that this would in many cases cause arrays to appear 
that contain heterogenous data.  In Julia, such arrays with abstract element
have severely lower performance.

The typical usage is
- Application developers assume inputs with consistent units.
- Application developers require any constants (acceleration of gravity, 
  gas constant...) as user input.
- Users convert all their input values as they define it `rho = 3 ← pound/foot^3`
- Users convert Muscade outputs before printing them out  
  `printf("stress [MPa] %f",stress → MPa)`

Inspect the file `espy.jl` for an overview of the available units.  

"""
    struct unit
        phydim :: Vector{Float64}
        ct     :: Float64
    end
    unit(phydim :: Vector{Float64})    = unit(phydim,1.)
    unit(ct     :: Float64 )           = unit([0.,0.,0.,0.,0.,0.,0.,0.],ct)

    Base.:(*)(a::unit  ,b::unit)       = unit( a.phydim+b.phydim,a.ct*b.ct)
    Base.:(*)(a::Number,b::unit)       = unit(          b.phydim,a   *b.ct)
    Base.:(*)(a::unit  ,b::Number)     = unit( a.phydim         ,a.ct*b   )
    Base.:(/)(a::unit  ,b::unit)       = unit( a.phydim-b.phydim,a.ct/b.ct)
    Base.:(/)(a::Number,b::unit)       = unit(         -b.phydim,a   /b.ct)
    Base.:(/)(a::unit  ,b::Number)     = unit( a.phydim         ,a.ct/b   )
    Base.:(^)(a::unit  ,b::Integer)    = unit( a.phydim*b       ,a.ct^b   )
    Base.:(^)(a::unit  ,b::Real)       = unit( a.phydim*b       ,a.ct^b   )
    Base.inv(a::unit)                  = unit(-a.phydim         ,1.  /a.ct)
    ←(a::Number,b::unit)               =                        a   *b.ct
    ←(a       ,b::unit)                =                        a  .*b.ct
    →(a::Number,b::unit)               =                        a   /b.ct
    →(a        ,b::unit)               =                        a  ./b.ct
    function string(a::unit)
        fundamental = ["m","kg","s","A","K","cd","mol","bit"]
        dsc =  a.ct == 1 ? "" : "$(a.ct) "
        for (i,d) in enumerate(a.phydim)
            if d==0.  continue end
            if d==1.  dsc = "$dsc$(fundamental[i])*"
            else      dsc = "$dsc$(fundamental[i])^$(d)*"
            end
        end
        if length(dsc)==0 dsc = "." end
        if dsc[end]=='*' dsc = dsc[1:end-1] end
        return dsc
    end
    show(io::IO,x::unit) = write(io,string(x))

    # Basic units
    const m        = metre    = unit([1.,0.,0.,0.,0.,0.,0.,0.])
    const kg       = kilogram = unit([0.,1.,0.,0.,0.,0.,0.,0.])
    const s        = second   = unit([0.,0.,1.,0.,0.,0.,0.,0.])
    const A        = Ampere   = unit([0.,0.,0.,1.,0.,0.,0.,0.])
    const K        = Kelvin   = unit([0.,0.,0.,0.,1.,0.,0.,0.])
    const Cd       = candela  = unit([0.,0.,0.,0.,0.,1.,0.,0.])
    const mol      = mole     = unit([0.,0.,0.,0.,0.,0.,1.,0.])
    const nat      = nit      = unit([0.,0.,0.,0.,0.,0.,0.,1.])

    # Prefixes
    const yocto    = unit(1e-24)
    const zepto    = unit(1e-21)
    const atto     = unit(1e-18)
    const femto    = unit(1e-15)
    const pico     = unit(1e-12)
    const nano     = unit(1e-9)
    const micro    = unit(1e-6)
    const milli    = unit(1e-3)
    const centi    = unit(1e-2)
    const deci     = unit(1e-1)
    const dimensionless = ena = unit(1e0)
    const deca     = unit(1e1)
    const hecto    = unit(1e2)
    const kilo     = unit(1e3)
    const mega     = unit(1e6)
    const giga     = unit(1e9)
    const tera     = unit(1e12)
    const peta     = unit(1e15)
    const exa      = unit(1e18)
    const zetta    = unit(1e21)
    const yotta    = unit(1e24)

    # Engineering
    const Å        = Angstrom = 1e-10metre
    const μm       = micrometre = micro*metre
    const mm       = milli*metre
    const cm       = centi*metre
    const dm       = deci*metre
    const km       = kilo*metre
    const are      = hecto*metre^2
    const ha       = hectare = hecto*are
    const l        = litre   = milli*metre^3
    const g        = gram    = milli*kilogram
    const Mg       = tonne = mega*gram
    const N        = Newton = m*kg/s^2
    const kN       = kilo*Newton
    const MN       = mega*Newton
    const GN       = giga*Newton
    const Pa       = Pascal = N/m^2
    const kPa      = kilo*Pascal
    const MPa      = mega*Pascal
    const GPa      = giga*Pascal
    const J        = Joule = N*m
    const W        = Watt = J/s
    const V        = Volt  = W/A #m^2*kg/s^3/A
    const mV       = milli*Volt
    const ε        = strain = dimensionless
    const με       = microstrain = micro
    const event    = dimensionless
    const Hz       = Hertz = 1/s
    const rad      = radian = dimensionless
    const period   = turn = 2π*dimensionless
    const deg      = degree = period/360
    const sr       = steradian = dimensionless
    const sphere   = 4π*sr
    const Coulomb  = C = A*s
    const Ohm      = Ω = V/A

    # Physical constants
    const c        = 299792458m/s
    const G        = 9.80665m/s^2
    const elementarycharge = e = 1.602176620898e-19C
    const Avogadro = Nₐ = 6.02214086e23/mol
    const Faraday  = F  = e*Nₐ
    const gasconst = R  = 8.314459848J/mol/K
    const Boltzman = kᵦ = R/Nₐ
    const Planck   = h  = 6.62601015e-35J*s
    const rPlanck  = ħ  = h/2π

    # data/probability
    const bit      = Shannon       = log( 2)nit
    const ban      = dit = Hartley = log(10)nit
    const octet    = byte = 8bit

    # Calendar
    const minute   = 60s
    const hour     = 60minute
    const day      = 24hour
    const week     = 7day
    const year     = 365.25day
    const month    = year/12

    # Prehistoric stuff. The imperials strike back. We can even do Donald Duck units!!!
    const eV       = elementarycharge*V
    const inch     = 25.4mm
    const foot     = 12inch
    const yard     = 3foot
    const fathom   = 2yard
    const furlong  = 110fathom
    const cable    = 120fathom
    const mils     = thou = milli*inch
    const nautical = 1852m
    const mile     = 1760yard
    const USacre   = 43560foot^2
    const mål      = 1000m^2
    const ounce    = 28.349523125g
    const pound    = 16ounce
    const shortton = 2kilo*pound
    const poundforce = pound*G
    const kip      = 1e3poundforce
    const cmil     = pi/4*mils^2
    const kcmil    = MCM = kilo*cmil
    const kgf      = kilogram*G
    const horsepower = 75kgf*m/s
    const knot     = nautical/hour
    const USpint   = 28.875*inch^3
    const USgallon = 8USpint
    const USfloz   = USpint/16
    const barrel   = 42USgallon
    const hogshead = 63USgallon
    const psi      = poundforce/inch^2
    const ksi      = kilo*psi
    const bar      = 1e5Pa
    const atm      = 101325Pa
    const mmHg     = 133.322387415Pa
    const torr     = (1/760)atm
    const BTU      = 1.0545e3J
    const calorie  = 4.184J
    const kgTNT    = 1e6calorie
    const alen     = 622.77cm
    const kWh      = 1e3W*hour
    const denier   = 1g/9000m