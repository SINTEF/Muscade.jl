#=
This is bugged
able_of_finite_diff_kernels is correct, table_of_finite_diff_kernels_transposed is apparently not, for 2nd derivatives

I need an algo that transposes the finite diff operator, and hashes a data structure that can be interogated using "finitediff"

=#

# table...[order+1][left or center or right][point in kernel]
const table_of_finite_diff_kernels            = [ [[(Δs=0,w=1.)] ],
                                                  [[(Δs=0,w=-1.),(Δs=1,w=1.)], [(Δs=-1,w=-1.),(Δs=0,w=1.)], [(Δs=-1,w=-1/2),(Δs=1,w=1/2)] ],
                                                  [[(Δs=0,w=1.),(Δs=1,w=-2.),(Δs=2,w=1.)], [(Δs=-2,w=1.),(Δs=-1,w=-2.),(Δs=0,w=1.)], [(Δs=-1,w=1.),(Δs=0,w=-2.),(Δs=1,w=1.)]]  ]

const table_of_finite_diff_kernels_transposed = [ [ [(Δs=0,w=1.)] ],
                                                  [[(Δs=0,w=-1.),(Δs=1,w=-1/2)],[(Δs=-1,w=1.),(Δs=1,w=-1/2)],[(Δs=-1,w=1/2),(Δs=1,w=-1.)],[(Δs=-1,w=1/2),(Δs=0,w=1.)],[(Δs=-1,w=1/2),(Δs=1,w=-1/2)]],
                                                  [[(Δs=0,w=1.),(Δs=1,w=1.)], [(Δs=-1,w=-2.),(Δs=0,w=-2.),(Δs= 1,w=1.)], [(Δs=-2,w=1.),(Δs=-1,w=1.),(Δs=0,w=-2.),(Δs= 1,w=1.)], [(Δs= 2,w=1.),(Δs= 1,w=1.),(Δs=0,w=-2.),(Δs=-1,w=1.)],
                                                   [(Δs= 1,w=-2.),(Δs=0,w=-2.),(Δs=-1,w=1.)], [(Δs=0,w=1.),(Δs=-1,w=1.)], [(Δs=-1,w=1.),(Δs=0,w=-2.),(Δs=1,w=1.)]                                            ] ] 
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
    #         x′[s+Δs] += x[s   ] * w/Δt  
    #     end
    #     for (Δs,w) ∈ finitediff(order,n,s)
    #         x′[s   ] += x[s+Δs] * w/Δt  
    #     end
    # end
    n<6 && Muscadeerror("Number of steps must be ≥6")
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

function FDsparsity(order,nstep)
    if     order == 0  nnz =   nstep
    elseif order == 1  nnz = 5*nstep-6
    elseif order == 2  nnz = 5*nstep-2 
    end
    istep = 𝕫1(undef,nnz)    
    jstep = 𝕫1(undef,nnz)
    if order ≥ 0
        istep[        1:  nstep  ] .= 1:nstep
        jstep[        1:  nstep  ] .= 1:nstep
    end
    if order ≥ 1
        istep[  nstep+1:2*nstep-2] .= 1:nstep-2
        jstep[  nstep+1:2*nstep-2] .= 3:nstep
        istep[2*nstep-1:3*nstep-4] .= 3:nstep
        jstep[2*nstep-1:3*nstep-4] .= 1:nstep-2
        istep[3*nstep-3:4*nstep-5] .= 1:nstep-1
        jstep[3*nstep-3:4*nstep-5] .= 2:nstep
        istep[4*nstep-4:5*nstep-6] .= 2:nstep
        jstep[4*nstep-4:5*nstep-6] .= 1:nstep-1
    end
    if order ≥ 2
        istep[5*nstep-5:5*nstep-2] .= [1,4,nstep-3,nstep]
        jstep[5*nstep-5:5*nstep-2] .= [4,1,nstep,nstep-3]
    end
    return istep,jstep         
end

