## FPL exact algorithm (2006)

#=
Based on the file:
http://ls11-www.cs.uni-dortmund.de/_media/rudolph/hypervolume/hv_python.zip
originally written in python by Simon Wessing (2010).

This implementation is for the variant 3 detailed in:
C. M. Fonseca, L. Paquete, and M. Lopez-Ibanez. An improved dimension-sweep
algorithm for the hypervolume indicator. In IEEE Congress on Evolutionary
Computation, pages 1157-1163, Vancouver, Canada, July 2006.
=#

mutable struct Node
    cargo::Vector
    next::Vector{Union{Nothing,Node}}
    prev::Vector{Union{Nothing,Node}}
    ignore::Int
    area::Vector{Float64}
    volume::Vector{Float64}
end

#=
function Base.show(io::IO, node::Node)
    print(io, "Node( ")
    print(io, "cargo = ", node.cargo)
    print(io, ", ")
    for n in node.next
        if isnothing(n)
            print(io, n, " -> ")
            continue
        end

        print(io, n.cargo, " -> ")
    end
    print(io, ")")

end
=#

function Node(numberLists, cargo = zeros(0);
        next   = fill(nothing, numberLists),
        prev   = fill(nothing, numberLists),
        ignore = 0,
        area   = zeros(numberLists),
        volume = zeros(numberLists),
    )

    Node(cargo, next, prev, ignore, area, volume)

end

function extend!(sentinel::Node, nodes, index)
    for node in nodes
        lastButOne = sentinel.prev[index]
        node.next[index] = sentinel
        node.prev[index] = lastButOne
        # set the last element as the new one
        sentinel.prev[index] = node
        lastButOne.next[index] = node
    end
end



"""
    Removes and returns 'node' from all lists in 1,2,..., 'index'.
"""
function remove!(node::Node, index, bounds)
    for i in 1:index-1
        predecessor = node.prev[i]
        successor = node.next[i]
        predecessor.next[i] = successor
        successor.prev[i] = predecessor

        if bounds[i] > node.cargo[i]
            bounds[i] = node.cargo[i]
        end
    end
    return node
end


"""
Inserts 'node' at the position it had in all lists in [1, 'index']
before it was removed. This method assumes that the next and previous
nodes of the node that is reinserted are in the list.

"""
function reinsert!(node::Node, index, bounds)
    for i in 1:index
        node.prev[i].next[i] = node
        node.next[i].prev[i] = node
        if bounds[i] > node.cargo[i]
            bounds[i] = node.cargo[i]
        end
    end

    node
end

"""
Sets up the list data structure needed for calculation.
"""
function preProcess(front, referencePoint)
    dimensions = length(referencePoint)
    sentinel = Node(dimensions)
    sentinel.next = fill(sentinel, dimensions)
    sentinel.prev = fill(sentinel, dimensions)
    nodes = [Node(dimensions, point) for point in front]

    for i in 1:dimensions
        # sortByDimension i
        sort!(nodes, by = node -> node.cargo[i])
        extend!(sentinel,nodes, i)
    end

    return sentinel

end

"""Recursive call to hypervolume calculation.

In contrast to the paper, the code assumes that the reference point
is [0, ..., 0]. This allows the avoidance of a few operations.
"""
function hvRecursive(sentinel, dimIndex, len, bounds)
    hvol = 0.0
    if len == 0
        return hvol
    elseif dimIndex == 1
        # special case: only one dimension
        # why using hypervolume at all?
        return -sentinel.next[1].cargo[1]

    elseif dimIndex == 2
        # special case: two dimensions, end recursion
        q = sentinel.next[2]
        h = q.cargo[1]
        p = q.next[2]

        while !(p === sentinel) #p is not sentinel
            pCargo = p.cargo
            hvol += h * (q.cargo[2] - pCargo[2])
            if pCargo[1] < h
                h = pCargo[1]
            end
            q = p
            p = q.next[2]
        end
        hvol += h * q.cargo[2]
        return hvol
    else
        p = sentinel
        q = p.prev[dimIndex]

        while !isempty(q.cargo) #q.cargo != None
            if q.ignore < dimIndex-1
                q.ignore = 0
            end
            q = q.prev[dimIndex]
        end

        q = p.prev[dimIndex]
        while len > 1 && (
                          q.cargo[dimIndex] > bounds[dimIndex] || q.prev[dimIndex].cargo[dimIndex] >= bounds[dimIndex])
            p = q
            remove!(p, dimIndex, bounds)
            q = p.prev[dimIndex]
            len -= 1
        end

        qArea = q.area
        qCargo = q.cargo
        qPrevDimIndex = q.prev[dimIndex]

        if len > 1
            hvol = qPrevDimIndex.volume[dimIndex]
            hvol += qPrevDimIndex.area[dimIndex] * ( qCargo[dimIndex] - qPrevDimIndex.cargo[dimIndex])
        else
            qArea[1] = 1
            qArea[2:dimIndex] = [qArea[i] * -qCargo[i] for i in 1:dimIndex-1]
        end

        q.volume[dimIndex] = hvol
        if q.ignore >= dimIndex-1
            qArea[dimIndex] = qPrevDimIndex.area[dimIndex]
        else

            qArea[dimIndex] = hvRecursive(sentinel, dimIndex - 1, len, bounds)
            if qArea[dimIndex] <= qPrevDimIndex.area[dimIndex]
                q.ignore = dimIndex-1
            end
        end

        while !(p === sentinel) #p is not sentinel
            pCargoDimIndex = p.cargo[dimIndex]
            hvol += q.area[dimIndex] * (pCargoDimIndex - q.cargo[dimIndex])
            bounds[dimIndex] = pCargoDimIndex
            reinsert!(p, dimIndex, bounds)
            len += 1
            q = p
            p = p.next[dimIndex]
            q.volume[dimIndex] = hvol

            if q.ignore >= dimIndex-1
                q.area[dimIndex] = q.prev[dimIndex].area[dimIndex]
            else
                q.area[dimIndex] = hvRecursive(sentinel, dimIndex - 1, len, bounds)
                if q.area[dimIndex] <= q.prev[dimIndex].area[dimIndex]
                    q.ignore = dimIndex-1
                end
            end

        end
        hvol -= q.area[dimIndex] * q.cargo[dimIndex]
        return hvol
    end
end

function FPL(front, referencePoint)
    dimensions = length(referencePoint)

    relevantPoints = [fx - referencePoint for fx in front]
    sentinel = preProcess(relevantPoints, referencePoint)

    bounds = fill(-Inf, dimensions)

    return hvRecursive(sentinel, dimensions, length(relevantPoints), bounds)
end