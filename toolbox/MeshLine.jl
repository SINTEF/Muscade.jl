"""
    MeshLine!(model,topNode, azimuth, eltype, xSection, segLength, nel)

Create the mesh of a multi-segment line in a Muscade model using beam or bar elements.
Nodes and elements are generated along a straight line, contained in the (:t1,:t2) plane, and oriented based on the given azimuth angle in radian, pointing from top to bottom.

# Arguments
- `model`: The Muscade model
- `topNode`: The top node of the mesh line.
- `azimuth`: The azimuth angle (in radians) defining the direction of the line.
- `eltype`: The type of element to use for the mesh, for example [`Bar3D`](@ref) or [`EulerBeam3D`](@ref).
- `xSection`: A vector of cross-section properties for each segment, for example [`BeamCrossSection`](@ref) or [`AxisymmetricBarCrossSection`](@ref)
- `segLength`: A vector of lengths for each segment.
- `nel`: A vector specifying the number of elements in each segment.

# Keyword Arguments
- `kwargsToElement...`: Additional keyword arguments passed to the element constructor (e.g., orientation parameters for beams).

# Returns
- `nodeList`: A vector of vectors containing the node IDs for each segment.
- `elementList`: A vector of element IDs for the entire mesh.
- `nodeCoord`: A vector of matrices containing the coordinates of nodes for each segment.

# Notes
Assumes at least 2 segments. 
"""
function MeshLine!(model, topNode, azimuth, eltype, xSection, segLength, nel; kwargsToElement...)
    # Function written in a very explicit style for easy understanding

    nseg            =   length(nel)
    accLength       =   [0; cumsum(segLength)]
    nnodes          =   nel .+ 1

    # Calculate node coordinates
    topNodeCoord    = coord([model.nod[topNode]])[1]
    bottomNodeCoord = topNodeCoord .+ [accLength[end] * cos(azimuth), accLength[end] * sin(azimuth), 0]
    nodeCoord       = [hcat(
            bottomNodeCoord[1] .- cos(azimuth) .* (accLength[seg].+((1:nnodes[seg]).- 1)/(nnodes[seg]-1)*segLength[seg]),
            bottomNodeCoord[2] .- sin(azimuth) .* (accLength[seg].+((1:nnodes[seg]).- 1)/(nnodes[seg]-1)*segLength[seg]),
            bottomNodeCoord[3] .+ zeros(Float64, nnodes[seg], 1)) for seg = 1:nseg];

    # Make room
    firstNode       =   Vector{Muscade.NodID}(undef, nseg)
    lastNode        =   Vector{Muscade.NodID}(undef, nseg)
    nodeList        =   Vector{Vector{Muscade.NodID}}(undef, nseg)
    elementList     =   Vector{Muscade.EleID}

    # Create Segment 1
    nodid           =   addnode!(model, nodeCoord[1])
    mesh            =   hcat(nodid[1:nnodes[1]-1], nodid[2:nnodes[1]])
    elementList     =   addelement!(model, eltype, mesh;                    mat=xSection[1], kwargsToElement...)
    firstNode[1]    =   nodid[1]
    lastNode[1]     =   nodid[size(nodid, 1)]
    nodeList[1]     =   nodid

    # Create intermediate segments (if they exist)
    if nseg > 2
        for segid ∈ 2:nseg-1
            local nodid         =   addnode!(model, nodeCoord[segid][2:end, :])
            firstNode[segid]    =   lastNode[segid-1]
            lastNode[segid]     =   nodid[size(nodid, 1)]
            local mesh          =   hcat(nodid[1:(nnodes[segid]-2)], nodid[2:(nnodes[segid]-1)])
            elementList         =   vcat(elementList, 
                addelement!(model, eltype, [firstNode[segid], nodid[1]];    mat=xSection[segid], kwargsToElement...), 
                addelement!(model, eltype, mesh;                            mat=xSection[segid], kwargsToElement...))
            nodeList[segid]     =   nodid
        end
    end

    # Create last segment
    nodid           =   addnode!(model, nodeCoord[nseg][2:end-1, :])
    firstNode[nseg] =   lastNode[nseg-1]
    lastNode[nseg]  =   topNode
    mesh            =   hcat(nodid[1:(nnodes[nseg]-3)], nodid[2:(nnodes[nseg]-2)])
    elementList     =   vcat(elementList, 
                addelement!(model, eltype, [firstNode[nseg], nodid[1]]; mat=xSection[nseg], kwargsToElement...),
                addelement!(model, eltype, mesh;                        mat=xSection[nseg], kwargsToElement...),
                addelement!(model, eltype, [nodid[end], topNode];       mat=xSection[nseg], kwargsToElement...))
    nodeList[nseg]  =   vcat(nodid, topNode)

    return nodeList, elementList, nodeCoord
end