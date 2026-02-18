# table...[order+1][left or center or right][point in kernel]
const table_of_finite_diff_kernels            = [ [[(Δs=0,w=1.)] ],
                                                  [[(Δs=0,w=-1.),(Δs=1,w=1.)], [(Δs=-1,w=-1.),(Δs=0,w=1.)], [(Δs=-1,w=-1/2),(Δs=1,w=1/2)] ],
                                                  [[(Δs=0,w=1.),(Δs=1,w=-2.),(Δs=2,w=1.)], [(Δs=-2,w=1.),(Δs=-1,w=-2.),(Δs=0,w=1.)], [(Δs=-1,w=1.),(Δs=0,w=-2.),(Δs=1,w=1.)]]  ]



function finitediff(order,n,s) 
    # INPUT
    # order - order of differentiation
    # n - length of series to differentiate
    # s - step index
    #
    # OUTPUT 
    # iterable collection of (Δs,w): index offset and weight
    #
    order>0 && n<6 && muscadeerror("Number of steps must be ≥6")
    p = table_of_finite_diff_kernels
    if     order == 0 p[1][1]
    elseif order == 1
        if     s==1   p[2][1]
        elseif s==n   p[2][2] 
        else          p[2][3]
        end
    elseif order == 2
        if     s==1   p[3][1]   
        elseif s==n   p[3][2] 
        else          p[3][3]  
        end
    end
end




