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
    function Base.string(a::unit)
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
    public const m        = metre    = unit([1.,0.,0.,0.,0.,0.,0.,0.])
    public const kg       = kilogram = unit([0.,1.,0.,0.,0.,0.,0.,0.])
    public const s        = second   = unit([0.,0.,1.,0.,0.,0.,0.,0.])
    public const A        = Ampere   = unit([0.,0.,0.,1.,0.,0.,0.,0.])
    public const K        = Kelvin   = unit([0.,0.,0.,0.,1.,0.,0.,0.])
    public const Cd       = candela  = unit([0.,0.,0.,0.,0.,1.,0.,0.])
    public const mol      = mole     = unit([0.,0.,0.,0.,0.,0.,1.,0.])
    public const nat      = nit      = unit([0.,0.,0.,0.,0.,0.,0.,1.]) # ... yes, but this *should* be a SI unit!

    # Prefixes
    public const yocto    = unit(1e-24)
    public const zepto    = unit(1e-21)
    public const atto     = unit(1e-18)
    public const femto    = unit(1e-15)
    public const pico     = unit(1e-12)
    public const nano     = unit(1e-9)
    public const micro    = unit(1e-6)
    public const milli    = unit(1e-3)
    public const centi    = unit(1e-2)
    public const deci     = unit(1e-1)
    public const dimensionless = ena = unit(1e0)
    public const deca     = unit(1e1)
    public const hecto    = unit(1e2)
    public const kilo     = unit(1e3)
    public const mega     = unit(1e6)
    public const giga     = unit(1e9)
    public const tera     = unit(1e12)
    public const peta     = unit(1e15)
    public const exa      = unit(1e18)
    public const zetta    = unit(1e21)
    public const yotta    = unit(1e24)

    # Engineering
    public const Å        = Angstrom = 1e-10metre
    public const μm       = micrometre = micro*metre
    public const mm       = milli*metre
    public const cm       = centi*metre
    public const dm       = deci*metre
    public const km       = kilo*metre
    public const are      = hecto*metre^2
    public const ha       = hectare = hecto*are
    public const l        = litre   = milli*metre^3
    public const g        = gram    = milli*kilogram
    public const Mg       = tonne = mega*gram
    public const N        = Newton = m*kg/s^2
    public const kN       = kilo*Newton
    public const MN       = mega*Newton
    public const GN       = giga*Newton
    public const Pa       = Pascal = N/m^2
    public const kPa      = kilo*Pascal
    public const MPa      = mega*Pascal
    public const GPa      = giga*Pascal
    public const J        = Joule = N*m
    public const W        = Watt = J/s
    public const V        = Volt  = W/A #m^2*kg/s^3/A
    public const mV       = milli*Volt
    public const ε        = strain = dimensionless
    public const με       = microstrain = micro
    public const event    = dimensionless
    public const Hz       = Hertz = 1/s
    public const rad      = radian = dimensionless
    public const period   = turn = 2π*dimensionless
    public const deg      = degree = period/360
    public const sr       = steradian = dimensionless
    public const sphere   = 4π*sr
    public const Coulomb  = C = A*s
    public const Ohm      = Ω = V/A

    # Physical constants
    public const c        = 299792458m/s
    public const G        = 9.80665m/s^2
    public const elementarycharge = e = 1.602176620898e-19C
    public const Avogadro = Nₐ = 6.02214086e23/mol
    public const Faraday  = F  = e*Nₐ
    public const gasconst = R  = 8.314459848J/mol/K
    public const Boltzman = kᵦ = R/Nₐ
    public const Planck   = h  = 6.62601015e-35J*s
    public const rPlanck  = ħ  = h/2π

    # data/probability
    public const bit      = Shannon       = log( 2)nit
    public const ban      = dit = Hartley = log(10)nit
    public const octet    = byte = 8bit

    # Imperial and other non SI
    public const minute   = 60s
    public const hour     = 60minute
    public const day      = 24hour
    public const week     = 7day
    public const year     = 365.25day
    public const month    = year/12
    public const kmh      = km/h

    public const eV       = elementarycharge*V
    public const inch     = 25.4mm
    public const foot     = 12inch
    public const yard     = 3foot
    public const fathom   = 2yard
    public const furlong  = 110fathom
    public const cable    = 120fathom
    public const mils     = thou = milli*inch
    public const nautical = 1852m
    public const mile     = 1760yard
    public const USacre   = 43560foot^2
    public const mål      = 1000m^2
    public const ounce    = 28.349523125g
    public const pound    = 16ounce
    public const shortton = 2kilo*pound
    public const poundforce = pound*G
    public const kip      = 1e3poundforce
    public const cmil     = pi/4*mils^2
    public const kcmil    = MCM = kilo*cmil
    public const kgf      = kilogram*G
    public const horsepower = 75kgf*m/s
    public const knot     = nautical/hour
    public const USpint   = 28.875*inch^3
    public const USgallon = 8USpint
    public const USfloz   = USpint/16
    public const barrel   = 42USgallon
    public const hogshead = 63USgallon
    public const psi      = poundforce/inch^2
    public const ksi      = kilo*psi
    public const bar      = 1e5Pa
    public const atm      = 101325Pa
    public const mmHg     = 133.322387415Pa
    public const torr     = (1/760)atm
    public const BTU      = 1.0545e3J
    public const calorie  = 4.184J
    public const kgTNT    = 1e6calorie
    public const alen     = 622.77cm
    public const kWh      = 1e3W*hour
    public const denier   = 1g/9000m