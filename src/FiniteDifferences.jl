
# table...[order+1][left or center or right][point in kernel]
const table_of_finite_diff_kernels            = [ [[(0,1.)] ],
                                                  [[(0,-1.),(1,1.)], [(-1,-1.),(0,1.)], [(-1,-1/2),(1,1/2)] ],
                                                  [[(0,1.),(1,-2.),(2,1.)], [(-2,1.),(-1,-2.),(0,1.)], [(-1,1.),(0,-2.),(1,1.)]]  ]

const table_of_finite_diff_kernels_transposed = [ [ [(0,1.)] ],
                                                  [[(0,-1.),(1,-1/2)],[(-1,1.),(1,-1/2)],[(-1,1/2),(1,-1.)],[(-1,1/2),(0,1.)],[(-1,1/2),(1,-1/2)]],
                                                  [[(0,1.),(1,1.)], [(-1,-2.),(0,-2.),( 1,1.)], [(-2,1.),(-1,1.),(0,-2.),( 1,1.)], [( 2,1.),( 1,1.),(0,-2.),(-1,1.)],
                                                   [( 1,-2.),(0,-2.),(-1,1.)], [(0,1.),(-1,1.)], [(-1,1.),(0,-2.),(1,1.)]                                            ] ] 
# npoints_transposed[order] NB: this is for a linear combination of derivatives up to `order`
function number_of_findiff_points(nstep,order)
   if     order == 1    nstep
   elseif order == 2  3*nstep-2
   elseif order == 3  3*nstep      
   end
end

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
    #         x′[s+Δs] += x[s   ] * w/Δt  
    #     end
    #     for (Δs,w) ∈ finitediff(order,n,s)
    #         x′[s   ] += x[s+Δs] * w/Δt  
    #     end
    # end
    p = table_of_finite_diff_kernels
    q = table_of_finite_diff_kernels_transposed
    if transposed
        if     order == 0 q[1][1]
        elseif order == 1
            if     s==1   q[2][1] 
            elseif s==2   q[2][2] 
            elseif s==n-1 q[2][3]
            elseif s==n   q[2][4]
            else          q[2][5] 
            end
        elseif order == 2 
            if     s==1   q[3][1]
            elseif s==2   q[3][2]   
            elseif s==3   q[3][3]
            elseif s==n-2 q[3][4]
            elseif s==n-1 q[3][5]
            elseif s==n   q[3][6]
            else          q[3][7]
            end
        end        
    else
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
end



