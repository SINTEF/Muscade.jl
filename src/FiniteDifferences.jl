

# table...[order+1][left or center or right][point in kernel]
const table_of_finite_diff_kernels            = [ [[(Δs=0,w=1.)]],
                                                  [[(Δs=0,w=-1/2),(Δs=2,w=1/2)], [(Δs=-2,w=-1/2),(Δs=0,w=1/2)], [(Δs=-1,w=-1/2),(Δs=1,w=1/2)] ],
                                                  [[(Δs=0,w=1.),(Δs=1,w=-2.),(Δs=2,w=1.)], [(Δs=-2,w=1.),(Δs=-1,w=-2.),(Δs=0,w=1.)], [(Δs=-1,w=1.),(Δs=0,w=-2.),(Δs=1,w=1.)]]  ]

const table_of_finite_diff_kernels_transposed = [ [ [(Δs=0,w=1.)] ],
                                                  [[(Δs=0,w=-.5),(Δs=1,w=-.5)], [(Δs=1,w=-.5)], [(Δs=-2,w=.5),(Δs=-1,w=.5),(Δs= 1,w=-.5)],
                                                  [(Δs= 2,w=-.5),(Δs=1,w=-.5),(Δs=-1,w=.5)], [(Δs=-1,w=.5)],[(Δs=0,w=.5),(Δs=-1,w=.5)],
                                                   [(Δs= -1,w=.5),(Δs= 1,w=-.5)]],
                                                  [[(Δs=0,w=1.),(Δs=1,w=1.)], [(Δs=-1,w=-2.),(Δs=0,w=-2.),(Δs= 1,w=1.)], [(Δs=-2,w=1.),(Δs=-1,w=1.),(Δs=0,w=-2.),(Δs= 1,w=1.)], 
                                                   [(Δs= 2,w=1.),(Δs= 1,w=1.),(Δs=0,w=-2.),(Δs=-1,w=1.)],[(Δs= 1,w=-2.),(Δs=0,w=-2.),(Δs=-1,w=1.)], [(Δs=0,w=1.),(Δs=-1,w=1.)], 
                                                   [(Δs=-1,w=1.),(Δs=0,w=-2.),(Δs=1,w=1.)]                                            ] ] 
# npoints_transposed[order] NB: this is for a linear combination of derivatives up to `order`



function finitediff(order,n,s;transposed=false) 
    # INPUT
    # order - order of differentiation
    # n - length of series to differentiate
    # s - step index
    # ;[transpose=false] 
    #
    # OUTPUT 
    # iterable collection of (Δs,w): index offset and weight
    #
    # transpose:
    # for s = 1:n
    #     for (Δs,w) ∈ finitediff(order,n,s,transpose=true)
    #         x′[s+Δs] += x[s   ] * w/Δt^order  
    #     end
    #     for (Δs,w) ∈ finitediff(order,n,s)
    #         x′[s   ] += x[s+Δs] * w/Δt^order  
    #     end
    # end
    n<6 && muscadeerror("Number of steps must be ≥6")
    p = table_of_finite_diff_kernels
    q = table_of_finite_diff_kernels_transposed
    if transposed
        if     order == 0 q[1][1]
        else
            if     s==1   q[order+1][1]
            elseif s==2   q[order+1][2]   
            elseif s==3   q[order+1][3]
            elseif s==n-2 q[order+1][4]
            elseif s==n-1 q[order+1][5]
            elseif s==n   q[order+1][6]
            else          q[order+1][7]
            end
        end        
    else
        if     order == 0 p[1][1]
        else
            if     s==1   p[order+1][1]   
            elseif s==n   p[order+1][2] 
            else          p[order+1][3]  
            end
        end
    end
end



