"""
    MeshLine(topNode, azimuth, eltype, xSection, segLength, nel)

Create the mesh of a multi-segment line, to initialize a finite element model.
Nodes and elements are generated along a straight line, contained in the (:t1,:t2) plane, and oriented based on the given azimuth angle.

# Arguments
- `topNode`: The top node of the mesh line.
- `azimuth`: The azimuth angle (in radians) defining the direction of the line.
- `eltype`: The type of element to use for the mesh, for example [`Bar3D`](@ref) or [`EulerBeam3D`](@ref).
- `xSection`: A vector of cross-section properties for each segment, for example [`BeamCrossSection`](@ref) or [`AxisymmetricBarCrossSection`](@ref)
- `segLength`: A vector of lengths for each segment.
- `nel`: A vector specifying the number of elements in each segment.

# Returns
- `nodeList`: A vector of vectors containing the node IDs for each segment.
- `elementList`: A vector of element IDs for the entire mesh.
- `nodeCoord`: A vector of matrices containing the coordinates of nodes for each segment.

# Notes
Assumes at least 2 segments.
"""
function MeshLine(topNode, azimuth, eltype, xSection, segLength, nel)
    # Need to pass element arguments (orient2 for EulerBeam3D for example)
    nseg = length(nel)
    accLength = [0; cumsum(segLength)]
    nnodes = nel .+ 1
    topNodeCoord = coord([topNode])[1]
    bottomNodeCoord = topNodeCoord .+ [accLength[end] * cos(azimuth), accLength[end] * sin(azimuth), 0]
    nodeCoord = [
        hcat(bottomNodeCoord[1] .- cos(azimuth) .* (accLength[seg] .+ ((1:nnodes[seg]) .- 1) / (nnodes[seg] - 1) * segLength[seg]),
            bottomNodeCoord[2] .- sin(azimuth) .* (accLength[seg] .+ ((1:nnodes[seg]) .- 1) / (nnodes[seg] - 1) * segLength[seg]),
            bottomNodeCoord[3] .+ zeros(Float64, nnodes[seg], 1)) for seg = 1:nseg
    ]

    # Lists with First and last node of each segment, etc.
    firstNode = Vector{Muscade.NodID}(undef, nseg)
    lastNode = Vector{Muscade.NodID}(undef, nseg)
    nodeList = Vector{Vector{Muscade.NodID}}(undef, nseg)
    elementList = Vector{Muscade.EleID}

    # Populate lists for Segment 1
    nodid = addnode!(model, nodeCoord[1])
    mesh = hcat(nodid[1:nnodes[1]-1], nodid[2:nnodes[1]])
    elementList = addelement!(model, eltype, mesh; mat=xSection[1])
    firstNode[1] = nodid[1]
    lastNode[1] = nodid[size(nodid, 1)]
    nodeList[1] = nodid

    # Populate list for the intermediate segments (if they exist)
    if nseg > 2
        for segid ∈ 2:nseg-1
            local nodid = addnode!(model, nodeCoord[segid][2:end, :])
            firstNode[segid] = lastNode[segid-1]
            lastNode[segid] = nodid[size(nodid, 1)]
            local mesh = hcat(nodid[1:(nnodes[segid]-2)], nodid[2:(nnodes[segid]-1)])
            elementList = vcat(elementList, addelement!(model, eltype, [firstNode[segid], nodid[1]]; mat=xSection[segid]))
            elementList = vcat(elementList, addelement!(model, eltype, mesh; mat=xSection[segid]))
            nodeList[segid] = nodid
        end
    end

    # Populate list for last segment
    nodid = addnode!(model, nodeCoord[nseg][2:end-1, :])
    firstNode[nseg] = lastNode[nseg-1]
    lastNode[nseg] = topNode.ID
    mesh = hcat(nodid[1:(nnodes[nseg]-3)], nodid[2:(nnodes[nseg]-2)])
    elementList = vcat(elementList, addelement!(model, eltype, [firstNode[nseg], nodid[1]]; mat=xSection[nseg]))
    elementList = vcat(elementList, addelement!(model, eltype, mesh; mat=xSection[nseg]))
    elementList = vcat(elementList, addelement!(model, eltype, [nodid[end], topNode.ID]; mat=xSection[nseg]))
    nodeList[nseg] = vcat(nodid, topNode.ID)

    return nodeList, elementList, nodeCoord
end