#module Performance

using Profile,ProfileView
using BenchmarkTools
using Muscade
using StaticArrays
using Muscade.Toolbox

N = 12
x = SVector{N,𝕣}(1:N)
X = revariate{2}(x)
i = SVector(1,2,3)
v = X[i]

apply = Muscade.apply{:chainrule}
#apply = Muscade.apply{:direct}

mission = :time
if mission == :eval
    out = apply(Toolbox.Rodrigues,v)
elseif mission == :time
    out =  apply(Toolbox.Rodrigues,v);
    @btime apply($Toolbox.Rodrigues,$v)
#    @btime Toolbox.Rodrigues($v)
elseif mission == :profile
    apply(Toolbox.Rodrigues,v)
    Profile.clear()
    Profile.@profile for i=1:250000
        local out = apply(Toolbox.Rodrigues,v)
    end
    ProfileView.view(fontsize=30);
    # After clicking on a bar in the flame diagram, you can type warntype_last() and see the result of 
    # code_warntype for the call represented by that bar.
end
;

#end # module


