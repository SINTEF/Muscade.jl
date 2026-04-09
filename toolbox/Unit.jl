    struct Unit
        phydim :: Vector{Float64}
        ct     :: Float64
    end
    Unit(phydim :: Vector{Float64})    = Unit(phydim,1.)
    Unit(ct     :: Float64 )           = Unit([0.,0.,0.,0.,0.,0.,0.,0.],ct)

    Base.:(*)(a::Unit  ,b::Unit)       = Unit( a.phydim+b.phydim,a.ct*b.ct)
    Base.:(*)(a::Number,b::Unit)       = Unit(          b.phydim,a   *b.ct)
    Base.:(*)(a::Unit  ,b::Number)     = Unit( a.phydim         ,a.ct*b   )
    Base.:(/)(a::Unit  ,b::Unit)       = Unit( a.phydim-b.phydim,a.ct/b.ct)
    Base.:(/)(a::Number,b::Unit)       = Unit(         -b.phydim,a   /b.ct)
    Base.:(/)(a::Unit  ,b::Number)     = Unit( a.phydim         ,a.ct/b   )
    Base.:(^)(a::Unit  ,b::Integer)    = Unit( a.phydim*b       ,a.ct^b   )
    Base.:(^)(a::Unit  ,b::Real)       = Unit( a.phydim*b       ,a.ct^b   )
    Base.inv(a::Unit)                  = Unit(-a.phydim         ,1.  /a.ct)
    ←(a::Number,b::Unit)               =                        a   *b.ct
    ←(a       ,b::Unit)                =                        a  .*b.ct
    →(a::Number,b::Unit)               =                        a   /b.ct
    →(a        ,b::Unit)               =                        a  ./b.ct
    function Base.string(a::Unit)
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
    show(io::IO,x::Unit) = write(io,string(x))

    # Basic units
    public m,metre,kg,kilogram,s,second,A,Ampere,K,Kelvin,Cd,candela,mol,mole,nat,nit
    const m        = metre    = Unit([1.,0.,0.,0.,0.,0.,0.,0.])
    const kg       = kilogram = Unit([0.,1.,0.,0.,0.,0.,0.,0.])
    const s        = second   = Unit([0.,0.,1.,0.,0.,0.,0.,0.])
    const A        = Ampere   = Unit([0.,0.,0.,1.,0.,0.,0.,0.])
    const K        = Kelvin   = Unit([0.,0.,0.,0.,1.,0.,0.,0.])
    const Cd       = candela  = Unit([0.,0.,0.,0.,0.,1.,0.,0.])
    const mol      = mole     = Unit([0.,0.,0.,0.,0.,0.,1.,0.])
    const nat      = nit      = Unit([0.,0.,0.,0.,0.,0.,0.,1.]) # ... yes, but this *should* be a SI Unit!

    # Prefixes
    public yocto,zepto,atto,femto,pico,nano,micromilli,centi,deci
    const yocto    = Unit(1e-24)
    const zepto    = Unit(1e-21)
    const atto     = Unit(1e-18)
    const femto    = Unit(1e-15)
    const pico     = Unit(1e-12)
    const nano     = Unit(1e-9)
    const micro    = Unit(1e-6)
    const milli    = Unit(1e-3)
    const centi    = Unit(1e-2)
    const deci     = Unit(1e-1)
    public dimensionless,ena,deca,hecto,kilo,mega,giga,tera,peta,exa,zetta,yotta
    const dimensionless = ena = Unit(1e0)
    const deca     = Unit(1e1)
    const hecto    = Unit(1e2)
    const kilo     = Unit(1e3)
    const mega     = Unit(1e6)
    const giga     = Unit(1e9)
    const tera     = Unit(1e12)
    const peta     = Unit(1e15)
    const exa      = Unit(1e18)
    const zetta    = Unit(1e21)
    const yotta    = Unit(1e24)

    # Engineering
    public Å,Angstrom,μm,micrometre,mm,millimetre,cm,centimetre,dm,decimetre,km,kilometre
    const Å        = Angstrom = 1e-10metre
    const μm       = micrometre = micro*metre
    const mm       = milli*metre
    const cm       = centi*metre
    const dm       = deci*metre
    const km       = kilo*metre
    public are,ha,hectare,l,litre,g,gram,Mg,tonne,N,Newton,kN,MN,GN
    const are      = hecto*metre^2
    const ha       = hectare = hecto*are
    const l        = litre   = milli*metre^3
    const g        = gram    = milli*kilogram
    const Mg       = tonne = mega*gram
    const N        = Newton = m*kg/s^2
    const kN       = kilo*Newton
    const MN       = mega*Newton
    const GN       = giga*Newton
    public Pa,Pascal,kPa,MPa,GPa,J,Joule,W,Watt,V,Volt,mV,ε,strain,με,microstrain
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
    public event,Hz,Hertz,rad,radian,period,turn,deg,degree,sr,steradian,sphere,Coulomb,C,Ohm,Ω
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
    public c,G,elementarycharge,e,Avogadro,Nₐ,Faraday,F,gasconst,R,Boltzman,kᵦ,Planck,h,rPlanck,ħ
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
    public bit,Shannon,ban,dit,Hartely,octet,byte
    const bit      = Shannon       = log( 2)nit
    const ban      = dit = Hartley = log(10)nit
    const octet    = byte = 8bit

    # Imperial and other non SI
    public minute,hour,day,week,year,month,kmh
    const minute   = 60s
    const hour     = 60minute
    const day      = 24hour
    const week     = 7day
    const year     = 365.25day
    const month    = year/12
    const kmh      = km/h

    public eV,inch,foot,yard,fathom,furlong,cable,mils,thou,nautical,mile,alen
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
    const alen     = 622.77cm
    public USacre,mål,ouce,pound,shortton,poundforce,kip,cmil,kcmil,MCM
    const USacre   = 43560foot^2
    const mål      = 1000m^2
    const ounce    = 28.349523125g
    const pound    = 16ounce
    const shortton = 2kilo*pound
    const poundforce = pound*G
    const kip      = 1e3poundforce
    const cmil     = pi/4*mils^2
    const kcmil    = MCM = kilo*cmil
    public kgf,horsepower,knot,USpint,USgallon,USfloz,barrel,hogshead
    const kgf      = kilogram*G
    const horsepower = 75kgf*m/s
    const knot     = nautical/hour
    const USpint   = 28.875*inch^3
    const USgallon = 8USpint
    const USfloz   = USpint/16
    const barrel   = 42USgallon
    const hogshead = 63USgallon
    public psi,ksi,bar,atm,mmHg,torr
    const psi      = poundforce/inch^2
    const ksi      = kilo*psi
    const bar      = 1e5Pa
    const atm      = 101325Pa
    const mmHg     = 133.322387415Pa
    const torr     = (1/760)atm
    public BTU,calorie,kgTNT,KWh,denier
    const BTU      = 1.0545e3J
    const calorie  = 4.184J
    const kgTNT    = 1e6calorie
    const kWh      = 1e3W*hour
    const denier   = 1g/9000m