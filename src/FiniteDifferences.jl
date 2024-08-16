
# if I want to compute x′, list of x values I need, and corresponding weights
function getFD(order,n,s,Δt) 
    if order == 0
        (i=𝕫[0],w=𝕣[1])
    elseif order == 1
        if     s==1   (i=𝕫[0,1],w=𝕣[-1,1]./Δt)
        elseif s==n-0 (i=𝕫[-1,0],w=𝕣[-1,1]./Δt) 
        else          (i=𝕫[-1,1],w=𝕣[-.5,.5]./Δt)
        end
    elseif order == 2
        if     s==1   (i=𝕫[0,1,2],w=𝕣[1,-2,1]./Δt^2)
        elseif s==n-0 (i=𝕫[-2,-1,0],w=𝕣[1,-2,1]./Δt^2)
        else          (i=𝕫[-1,0,1],w=𝕣[1,-2,1]./Δt^2)
        end
    end
end
# I have a value of x, where do I add it in in x′ and with what ceofficients?
# Usefull if a value of x is available, but should no be stored
function gettransposedFD(order,n,s,Δt) 
    if order == 0
        (i=𝕫[0],w=𝕣[1])
    elseif order == 1
        if     s==1   (i=𝕫[0,1],w=𝕣[-1,-.5]./Δt)
        elseif s==2   (i=𝕫[-1,1],w=𝕣[1,-.5]./Δt)
        elseif s==n-1 (i=𝕫[-1,1],w=𝕣[.5,-1]./Δt) 
        elseif s==n-0 (i=𝕫[-1,0],w=𝕣[.5,1]./Δt) 
        else          (i=𝕫[-1,1],w=𝕣[.5,-.5]./Δt)
        end
    elseif order == 2
        if     s==1   (i=𝕫[0,1],w=𝕣[1,1]./Δt^2)
        elseif s==2   (i=𝕫[-1,0,1],w=𝕣[-2,-2,1]./Δt^2)
        elseif s==3   (i=𝕫[-2,-1,0,1],w=𝕣[1,1,-2,1]./Δt^2)
        elseif s==n-2 (i=𝕫[-1,0,1,2],w=𝕣[1,-2,1,1]./Δt^2)  
        elseif s==n-1 (i=𝕫[-1,0,1],w=𝕣[1,-2,-2]./Δt^2)  
        elseif s==n-0 (i=𝕫[-1,0],w=𝕣[1,1]./Δt^2)
        else          (i=𝕫[-1,0,1],w=𝕣[1,-2,1]./Δt^2)
        end
    end
end

