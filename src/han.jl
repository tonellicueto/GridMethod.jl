#  STRUCTURES

mutable struct Rectangle
    x_min::Float64
    x_max::Float64
    y_min::Float64
    y_max::Float64
    segs::Array{MySegment,1}
    Rectangle(x_min::Float64, x_max::Float64, y_min::Float64, y_max::Float64) = new(x_min, x_max, y_min, y_max)
end

mutable struct KDData
    reflex_vertices::Array{Int64,1}
    contour::Array{Vertex, 1}
    on::Array{Int64, 1}
    intH::Array{MySegment, 1}
    intV::Array{MySegment, 1}
    rectangles::Array{Rectangle,1}
    minBB::Rectangle
    KDData() = new()
end

mutable struct KDNode{N, T}
    boundary::HyperRectangle{N, T}
    data::KDData
    children::Union{Array{KDNode{N,T}, 1}, Void}
    parent::Union{KDNode{N, T}, Void}
    split_dir::Int64
    split_val::Float64
    visited::Bool
end

function KDNode(origin::SVector{N, T}, widths::SVector{N, T}, data::KDData=nothing) where {N, T}
    return KDNode{N,T}(HyperRectangle(origin, widths), data, nothing, nothing, 0, 0.0, false)
end

mutable struct Face
    contour::Array{Tuple{Float64, Float64}, 1}
    coord::Int64
    dir::Int64
    aff::Float64
    root::KDNode
    Face(contour::Array{Tuple{Float64, Float64}, 1}, coord::Int64, dir::Int64, aff::Float64 )= new(contour, coord, dir,aff)
end


function strictly_intersecting(t::Array{MySegment,1}, bd::HyperRectangle)
    output = MySegment[]
    for i = 1:length(t)
    s = t[i]
    a = bd.origin[1+s.orient]
    b = bd.origin[2-s.orient]
    if (b< s.aff && s.aff < b + bd.widths[2-s.orient]) &&
        min(s.v1, s.v2) < a + bd.widths[1+s.orient] && #abs(min(s.v1, s.v2) - a - bd.widths[1+s.orient] > 0 ) ) &&
        a < max(s.v1, s.v2) # || abs(a - max(s.v1, s.v2)) > 0 )
        #(a <= s.v1 && s.v1 <= a + cell.boundary.widths[1] || a <= s.v2 && s.v2 <= a + cell.boundary.widths[1]) || (b == s.aff &&
            push!(output, s)
    end
    end
    return output
end

function intersecting(t::Array{MySegment,1}, bd::Rectangle)
    output = MySegment[]
    for i = 1:length(t)
    s = t[i]
    origin = (bd.x_min, bd.y_min)
    widths = (bd.x_max - bd.x_min, bd.y_max - bd.y_min)
    a = origin[1+s.orient]
    b = origin[2-s.orient]
    if (b<= s.aff && s.aff <= b + widths[2-s.orient]) && min(s.v1, s.v2) < a + widths[1+s.orient] && a < max(s.v1, s.v2)
            push!(output, s)
    end
    end
    return output
end


#------------------------   KD TREE -------------------------------

@inline isleaf(node::KDNode) = node.children === nothing
@inline children(node::KDNode) = node.children
@inline parent(node::KDNode) = node.parent
@inline center(node::KDNode) = center(node.boundary)
@inline vertices(node::KDNode) = vertices(node.boundary)


function inside(x::Float64,y::Float64, bd::HyperRectangle)
    return bd.origin[1]<=x && x<=bd.origin[1]+bd.widths[1] && bd.origin[2]<=y && y<=bd.origin[2]+bd.widths[2]
end

function strictly_inside(x::Float64,y::Float64, bd::HyperRectangle)
    return bd.origin[1]<x && x<bd.origin[1]+bd.widths[1] && bd.origin[2]<y && y<bd.origin[2]+bd.widths[2]
end

function inside(x::Float64,y::Float64, bd::Rectangle)
    return bd.x_min<=x && x<=bd.x_max && bd.y_min<=y && y<=bd.y_max
end

function splitNode!(node::KDNode, v::Array{Vertex, 1})
    split_coord = node.split_dir
    node.split_val = mymedian([v[p].coords[split_coord+1] for p in node.data.reflex_vertices])
    split_val = node.split_val

    c = Array{KDNode,1}(2)
    if split_coord == 0 # split at x
        c[1] = KDNode(node.boundary.origin, SVector(split_val-node.boundary.origin[1], node.boundary.widths[2]), KDData())
        c[2] = KDNode(SVector(split_val, node.boundary.origin[2]), SVector(node.boundary.widths[1]-c[1].boundary.widths[1], node.boundary.widths[2]), KDData())
    else #split at y
        c[1] = KDNode(node.boundary.origin, SVector(node.boundary.widths[1],split_val-node.boundary.origin[2]), KDData())
        c[2] = KDNode(SVector(node.boundary.origin[1],split_val), SVector(node.boundary.widths[1],node.boundary.widths[2]-c[1].boundary.widths[2]), KDData())
    end
    c[1].split_dir = (split_coord+1)%2
    c[2].split_dir = (split_coord+1)%2
    c[1].data.reflex_vertices = Int64[]
    c[2].data.reflex_vertices = Int64[]
    c[1].data.rectangles = Rectangle[]
    c[2].data.rectangles = Rectangle[]
    node.children = [c[1], c[2]]
    for child in node.children
        child.parent = node
    end
    #----refine data in each child
    for j = 1:length(node.data.reflex_vertices)
        n = node.data.reflex_vertices[j]
        if v[n].coords[split_coord+1] < split_val
            push!(c[1].data.reflex_vertices, node.data.reflex_vertices[j])
        elseif v[n].coords[split_coord+1] > split_val
            push!(c[2].data.reflex_vertices, node.data.reflex_vertices[j])
        end
    end
end

function allnodes(node::KDNode)
    Channel() do c
        queue = [node]
        while !isempty(queue)
            current = pop!(queue)
            put!(c, current)
            if !isleaf(current)
                append!(queue, children(current))
            end
        end
    end
end

function allleaves(node::KDNode)
    Channel() do c
        for child in allnodes(node)
            if isleaf(child)
                put!(c, child)
            end
        end
    end
end

function mymedian(n)
	s = sort(n)
	len = length(n)
	return  s[floor(Int, len / 2) + 1]
end

function get_rect!(node::KDNode, E_h::Vector{MySegment}, E_v::Vector{MySegment})
    node.data.rectangles = Rectangle[]
    node.data.intH = MySegment[]
    node.data.intV = MySegment[]
    node.data.intH = append!(node.data.intH, strictly_intersecting(E_h, node.boundary))
    node.data.intV = append!(node.data.intV, strictly_intersecting(E_v, node.boundary))
    inter = MySegment[]
    indx = Int64[]
    SortEdges(node.data.intV)
    SortEdges(node.data.intH)
    if length(node.data.intV) + length(node.data.intH) == 0
        push!(node.data.rectangles, Rectangle(node.boundary.origin[1],node.boundary.origin[1] + node.boundary.widths[1], node.boundary.origin[2],node.boundary.origin[2] + node.boundary.widths[2] ))
        return
    end
    if node.parent == nothing
        push!(node.data.rectangles, Rectangle(node.data.intV[1].aff,node.data.intV[2].aff, node.data.intH[1].aff,node.data.intH[2].aff ))
        return
    end
        i = 1
        while i <= length(node.data.intH)
            s = node.data.intH[i]
            for j = 1:length(node.data.intV)
                 t= node.data.intV[j]
                if t.aff == s.v1 || t.aff == s.v2
                    push!(inter, t)
                    push!(indx, j)
                end
            end
            if length(inter) == 2
                deleteat!(node.data.intV, indx[1])
                deleteat!(node.data.intV, indx[2]-1)
                pe = 0.0
                if inter[1].v1 == s.aff
                    pe = inter[1].v2
                else
                    pe = inter[1].v1
                end
                if pe > s.aff
                    push!(node.data.rectangles, Rectangle(min(inter[1].aff, inter[2].aff), max(inter[1].aff, inter[2].aff),s.aff,node.boundary.origin[2] + node.boundary.widths[2] ))
                else
                    push!(node.data.rectangles, Rectangle(min(inter[1].aff, inter[2].aff), max(inter[1].aff, inter[2].aff), node.boundary.origin[2], s.aff))
                end
                deleteat!(node.data.intH, i)
                i = i-1
            end
            inter = MySegment[]
            indx = Int64[]
            i+=1
        end

        i = 1
        while i <= length(node.data.intV)
            s = node.data.intV[i]
            for j = 1:length(node.data.intH)
                 t= node.data.intH[j]
                if t.aff == s.v1 || t.aff == s.v2
                    push!(inter, t)
                    push!(indx, j)
                end
            end
            if length(inter) == 2
                deleteat!(node.data.intH, indx[1])
                deleteat!(node.data.intH, indx[2]-1)
                pe = 0.0
                if inter[1].v1 == s.aff
                    pe = inter[1].v2
                else
                    pe = inter[1].v1
                end
                if pe > s.aff
                    push!(node.data.rectangles, Rectangle(s.aff,node.boundary.origin[1] + node.boundary.widths[1],min(inter[1].aff, inter[2].aff), max(inter[1].aff, inter[2].aff)))
                else
                    push!(node.data.rectangles, Rectangle(node.boundary.origin[1], s.aff, min(inter[1].aff, inter[2].aff), max(inter[1].aff, inter[2].aff)))
                end
                deleteat!(node.data.intV, i)
                i = i-1
            end
            inter = MySegment[]
            indx = Int64[]
            i+=1
        end

        i = 1
        while i <= length(node.data.intH)
            s = node.data.intH[i]
            for j = 1:length(node.data.intV)
                 t= node.data.intV[j]
                if t.aff == s.v1 || t.aff == s.v2
                    push!(inter, t)
                    push!(indx, j)
                end
            end
            if length(inter) == 1
                deleteat!(node.data.intV, indx[1])
                pe = 0.0
                pe2 = 0.0
                if inter[1].v1 == s.aff
                    pe = inter[1].v2
                else
                    pe = inter[1].v1
                end
                if s.v1 == inter[1].aff
                    pe2 = s.v2
                else
                    pe2 = s.v1
                end
                if pe > s.aff && pe2 > inter[1].aff
                    push!(node.data.rectangles, Rectangle(inter[1].aff, node.boundary.origin[1] + node.boundary.widths[1]  ,s.aff,node.boundary.origin[2] + node.boundary.widths[2] ))
                elseif pe > s.aff && pe2 < inter[1].aff
                    push!(node.data.rectangles, Rectangle(node.boundary.origin[1], inter[1].aff, s.aff,node.boundary.origin[2] + node.boundary.widths[2] ))
                elseif pe < s.aff && pe2 > inter[1].aff
                    push!(node.data.rectangles, Rectangle(inter[1].aff, node.boundary.origin[1] + node.boundary.widths[1], node.boundary.origin[2], s.aff))
                elseif pe < s.aff && pe2 < inter[1].aff
                    push!(node.data.rectangles, Rectangle(node.boundary.origin[1], inter[1].aff, node.boundary.origin[2], s.aff))
                end
                deleteat!(node.data.intH, i)
                i = i-1
            end
            inter = MySegment[]
            indx = Int64[]
            i+=1
        end

        i = 1
        while i <= length(node.data.intH)
            s = node.data.intH[i]
            if i == 1 && s.v2 < s.v1
                push!(node.data.rectangles, Rectangle(node.boundary.origin[1],node.boundary.origin[1] + node.boundary.widths[1], node.boundary.origin[2],s.aff ))
                deleteat!(node.data.intH, i)
                i = i-1
            elseif i ==1 && s.v2 > s.v1
                if length(node.data.intH)>1
                    t = node.data.intH[2]
                    push!(node.data.rectangles, Rectangle(node.boundary.origin[1],node.boundary.origin[1] + node.boundary.widths[1], s.aff,t.aff ))
                    deleteat!(node.data.intH, i)
                    deleteat!(node.data.intH, i)
                    i = i-1
                else
                    push!(node.data.rectangles, Rectangle(node.boundary.origin[1],node.boundary.origin[1] + node.boundary.widths[1], s.aff, node.boundary.origin[2] + node.boundary.widths[2]))
                    deleteat!(node.data.intH, i)
                    i = i-1
                end
             end
             i+=1
        end

        i = 1
        while i <= length(node.data.intV)
            s = node.data.intV[i]
            if i == 1 && s.v1 < s.v2
                push!(node.data.rectangles, Rectangle(node.boundary.origin[1],s.aff, node.boundary.origin[2] , node.boundary.origin[2]+ node.boundary.widths[2]))
                deleteat!(node.data.intV, i)
                i = i-1
            elseif i ==1 && s.v1 > s.v2
                if length(node.data.intV)>1
                    t = node.data.intV[2]
                    push!(node.data.rectangles, Rectangle( s.aff,t.aff , node.boundary.origin[2] , node.boundary.origin[2]+ node.boundary.widths[2]))
                    deleteat!(node.data.intV, i)
                    deleteat!(node.data.intV, i)
                    i = i-1
                else
                    push!(node.data.rectangles, Rectangle(s.aff, node.boundary.origin[1] + node.boundary.widths[1], node.boundary.origin[2] , node.boundary.origin[2]+ node.boundary.widths[2]))
                    deleteat!(node.data.intV, i)
                    i = i-1
                end
             end
             i+=1
        end
end

function KDTree(root::KDNode, v::Array{Vertex, 1})
    queue = [root]
    #dir = 0
    while !isempty(queue)
        node = pop!(queue)
        if length(node.data.reflex_vertices) > 0
            splitNode!(node, v)
            append!(queue, children(node))
            #dir += 1
        end
    end
end

function ccw(p1::Tuple{Float64, Float64}, p2::Tuple{Float64, Float64}, p3::Tuple{Float64, Float64})
    res =  (p2[1] - p1[1])*(p3[2] - p1[2]) - (p2[2] - p1[2])*(p3[1] - p1[1])
    if res >= 0
        return true
    else
        return false
    end
end

function BV(R1::Rectangle, R2::Rectangle)
    return Rectangle(min(R1.x_min, R2.x_min), max(R1.x_max, R2.x_max), min(R1.y_min, R2.y_min), max(R1.y_max,R2.y_max))
end


function cover(R1::Rectangle, R2::Rectangle)
    # returns true if R1 covers R2
    #one rectangle is on the left side of other
        return R2.x_min >= R1.x_min && R2.y_min >= R1.y_min && R2.x_max <= R1.x_max && R2.y_max <= R1.y_max
end

function overlap(R1::Rectangle, R2::Rectangle)
    #one rectangle is on the left side of other
        if (R1.x_min > R2.x_max || R2.x_min > R1.x_max)
            return false
        end
    #one rectangle is above other
        if (R1.y_max < R2.y_min || R2.y_max < R1.y_min)
            return false
        end

        return true
end

function overlap(R1::Rectangle, bd::HyperRectangle)
    #one rectangle is on the left side of other
    xmin = bd.origin[1]
    xmax = bd.origin[1] + bd.widths[1]
    ymin = bd.origin[2]
    ymax = bd.origin[2] + bd.widths[2]
        if (R1.x_min > xmax || xmin > R1.x_max)
            return false
        end
    #one rectangle is above other
        if (R1.y_max < ymin || ymax < R1.y_min)
            return false
        end

        return true
end

#returns true if rectangle Q intesects the collection of rectangles of the decomposition
function rectangleIntersects(root::KDNode, Q::Rectangle)
    flag = false
    st = Stack(KDNode)
    push!(st,root)
    while !isempty(st)
        node = pop!(st)
        if cover(Q, node.data.minBB)
            flag = true
            return true
        end

        if overlap(Q, node.data.minBB)
            if isleaf(node)
                for r in node.data.rectangles
                    if overlap(Q,r)
                        flag = true
                        return true
                    end
                end
            else
                if overlap(Q, children(node)[1].data.minBB)
                    push!(st, children(node)[1])
                end

                if overlap(Q, children(node)[2].data.minBB)
                    push!(st, children(node)[2])
                end
            end
        else #not overlap with the minBB
            return false
        end
    end
    if flag != true
        return false
    end
end

#----------------needed for the 2d algorithm------

function pointIsOut(root::KDNode, p::Tuple{Float64,Float64})
    #true if point is outside, false otherwise
    st = Stack(KDNode)
    push!(st,root)
    while !isempty(st)
        node = pop!(st)
        if !inside(p[1], p[2], node.data.minBB)
            if isempty(st)
                return true
            end
        elseif isleaf(node)
            flag = true
            for r in node.data.rectangles
                if inside(p[1], p[2], r)
                    flag = false
                    break
                end
            end
            return flag
        else
            if p[node.split_dir+1] < node.split_val
                push!(st, children(node)[1])
            elseif p[node.split_dir+1] > node.split_val
                push!(st, children(node)[2])
            else
                push!(st, children(node)[1])
                push!(st, children(node)[2])
            end
        end
    end
    return false
end

function findIntersectingRects(root::KDNode, Q::Rectangle)
    output = Rectangle[]
    st = Stack(KDNode)
    push!(st,root)
    while !isempty(st)
        node = pop!(st)
        if cover(Q, node.data.minBB)
            for r in allleaves(node)
                append!(output, r.data.rectangles)
            end
        elseif overlap(Q, node.data.minBB)
            if isleaf(node)
                for r in node.data.rectangles
                    if overlap(Q,r)
                        push!(output, r)
                    end
                end
            else
                if overlap(Q, children(node)[1].boundary)
                    push!(st, children(node)[1])
                end

                if overlap(Q, children(node)[2].boundary)
                    push!(st, children(node)[2])
                end
            end
        end
    end
    return output
end


function set_rects!(root::KDNode, E_h::Vector{MySegment}, E_v::Vector{MySegment})
    i = 1
 for node in allleaves(root)
     #println(i)
        get_rect!(node, E_h, E_v)
        for rect in node.data.rectangles
        rect.segs = intersecting(E_v, rect)
        append!(rect.segs, intersecting(E_h, rect))
        end
        i+=1
  end # end of for loop
end

# initialize the search tree with the min bounding volumes
function searchTree!(root::KDNode)
    s = Stack(KDNode)
    push!(s,root)
    while !isempty(s)
        n = top(s)
        if !isleaf(n) && n.visited == false
            push!(s, children(n)[1])
            push!(s, children(n)[2])
            n.visited = true
        elseif !isleaf(n) && n.visited == true
            n.data.minBB = deepcopy(n.data.rectangles[1])
            for i = 2:length(n.data.rectangles)
                n.data.minBB = BV(n.data.minBB, n.data.rectangles[i] )
            end
            if n!= root
                push!(parent(n).data.rectangles, n.data.minBB)
            end
            pop!(s)
        elseif isleaf(n) && n.parent != nothing
            n.data.minBB = deepcopy(n.data.rectangles[1])
            for i = 2:length(n.data.rectangles)
                n.data.minBB = BV(n.data.minBB, n.data.rectangles[i] )
            end
            pop!(s)
            push!(parent(n).data.rectangles, n.data.minBB)
        elseif isleaf(n) && n.parent == nothing
            n.data.minBB = deepcopy(n.data.rectangles[1])
            pop!(s)
        end
    end
end

function LinfDist(p::Tuple{Float64,Float64}, s::MySegment)
    if p[s.orient+1] < min(s.v1, s.v2)
        return max( abs(p[s.orient+1] - min(s.v1, s.v2)), abs(p[2-s.orient]-s.aff))
    elseif p[s.orient+1] > max(s.v1, s.v2)
        return max( abs(p[s.orient+1] - max(s.v1, s.v2)), abs(p[2-s.orient]-s.aff))
    else
        return abs(p[2-s.orient]-s.aff)
    end
end

function estimateDist(root::KDNode, p::Tuple{Float64,Float64})
    rects = Rectangle[]
    dist = Inf
    st = Stack(KDNode)
    push!(st,root)
    while !isempty(st)
        node = pop!(st)
        if !inside(p[1], p[2], node.data.minBB)
            return dist
        end
        if isleaf(node)
            for r in node.data.rectangles
                if inside(p[1], p[2], r)
                    push!(rects, r)
                    flag = true
                end
            end
        else
            if p[node.split_dir+1] < node.split_val
                push!(st, children(node)[1])
            elseif p[node.split_dir+1] > node.split_val
                push!(st, children(node)[2])
            else
                push!(st, children(node)[1])
                push!(st, children(node)[2])
            end
        end
    end
    if length(rects) == 0
        return dist # dist is equal to Inf by now
    else
        for r in rects
            for s in r.segs
                if OrientedZone(s, p) && dist_from_seg(s,p) < dist # or LinfDist(p,s) < dist
                    dist = LinfDist(p,s)
                end
            end
        end
    end
    return dist
end

function decomposition!(F::Vector{Array{Tuple{Float64,Float64}}}, E::Vector{Vector{MySegment}})
    c= 0
    for e in E
        c += length(e)
    end
    v = Array{Vertex, 1}(c)
    E_v = Vector{MySegment}(Int(c/2))
    E_h = Vector{MySegment}(Int(c/2))
    c = 0
    k = length(E)
    for i=1:k
        v[c+1] = Vertex(F[i][1])
        v[c+1].indx = c+1
        for j=2:length(F[i])
            v[c+j] = Vertex(F[i][j])
            v[c+j].indx = c+j
            v[c+j-1].next = v[c+j]
            v[c+j].prev = v[c+j-1]
        end
        v[c+ length(F[i])].next = v[c+1]
        v[c+1].prev = v[c + length(F[i])]
        c += length(F[i])
    end

    c1 = 1
    c2 = 1
    for e in E
        for s in e
            if s.orient ==1
                E_v[c1] = s
                c1+=1
            else
                E_h[c2] = s
                c2 +=1
            end
        end
    end
    SortEdges(E_v)
    SortEdges(E_h)

    for i =1:length(v)
        if !ccw(v[i].prev.coords, v[i].coords, v[i].next.coords)# angle is concave
            v[i].reflex = true
        end
    end

    #decomposition
    d = KDData()
    d.reflex_vertices = Tuple{Float64, Float64, Int64}[]
    for i = 1:length(v)
        if v[i].reflex == true
            push!(d.reflex_vertices, i)
        end
    end

    d.intH = E_h
    d.intV = E_v
    d.contour = v
    d.rectangles = Rectangle[]
    d.on = Int64[]

 # define bounding box
    x_min = Inf
    y_min = Inf
    x_max = -Inf
    y_max = -Inf

    for f in E
        for s in f
            if s.orient && s.aff < x_min
                x_min = s.aff
            end
            if s.orient && s.aff > x_max
                x_max = s.aff
            end
            if !s.orient && s.aff < y_min
                y_min = s.aff
            end
            if !s.orient && s.aff > y_max
                y_max = s.aff
            end
        end
    end
    #x_min = E_v[1].aff
    #x_max = E_v[length(E_h)].aff
    #y_min = E_h[1].aff
    #y_max = E_h[length(E_h)].aff
    l = max(abs(x_max- x_min), abs(y_max - y_min) )

    root2 = KDNode( SVector(x_min-0.5, y_min-0.5), SVector(l+2.,l+2.), d)

    root2.split_dir = 0

    KDTree(root2,v)
    set_rects!(root2, E_h,E_v)
    searchTree!(root2)

    root2
end

#----- for the 3d algorithm
function decompose!(f::Face)
	E = MySegment[]
    E = getEdges(f.contour, 0)
    f.root = decomposition!(f.contour,E)
end

function decomposition!(F::Array{Tuple{Float64,Float64}}, E::Vector{MySegment})
    v = Array{Vertex}(length(F))
    E_v = Array{MySegment}(Int(length(F)/2))
    E_h = Array{MySegment}(Int(length(F)/2))
        v[1] = Vertex(F[1])
        v[1].indx = 1
        for j=2:length(F)
            v[j] = Vertex(F[j])
            v[j].indx = j
            v[j-1].next = v[j]
            v[j].prev = v[j-1]
        end
        v[length(F)].next = v[1]
        v[1].prev = v[length(F)]
    c1 = 1
    c2 = 1
    for s in E
            if s.orient ==1
                E_v[c1] = s
                c1+=1
            else
                E_h[c2] = s
                c2 +=1
            end
    end
    SortEdges(E_v)
    SortEdges(E_h)

    for i =1:length(v)
        if !ccw(v[i].prev.coords, v[i].coords, v[i].next.coords)# angle is concave
            v[i].reflex = true
        end
    end

    #decomposition
    d = KDData()
    d.reflex_vertices = Tuple{Float64, Float64, Int64}[]
    for i = 1:length(v)
        if v[i].reflex == true
            push!(d.reflex_vertices, i)
        end
    end

    d.intH = E_h
    d.intV = E_v
    d.contour = v
    d.rectangles = Rectangle[]
    d.on = Int64[]

 # define bounding box
    x_min = Inf
    y_min = Inf
    x_max = -Inf
    y_max = -Inf

    for s in E
            if s.orient && s.aff < x_min
                x_min = s.aff
            end
            if s.orient && s.aff > x_max
                x_max = s.aff
            end
            if !s.orient && s.aff < y_min
                y_min = s.aff
            end
            if !s.orient && s.aff > y_max
                y_max = s.aff
            end
    end
    #x_min = E_v[1].aff
    #x_max = E_v[length(E_h)].aff
    #y_min = E_h[1].aff
    #y_max = E_h[length(E_h)].aff
    l = max(abs(x_max- x_min), abs(y_max - y_min) )

    root2 = KDNode( SVector(x_min-0.5, y_min-0.5), SVector(l+2.,l+2.), d)

    root2.split_dir = 0

    KDTree(root2,v)
    set_rects!(root2, E_h,E_v)
    searchTree!(root2)

    root2
end

function SemiAlgebraicTypes.mesh(root::KDNode)
           m = SemiAlgebraicTypes.mesh(Float64)
           for leaf in allleaves(root)
               v = hcat(collect(vertices(leaf.boundary))...)
               if size(v,1) == 2
                   v = vcat(v,fill(0.,1,size(v,2)))
               end
               n = nbv(m)
               for j in 1:size(v,2)
                   push_vertex!(m, v[1:3,j])
               end
               push_edge!(m, [n+1,n+2])
               push_edge!(m, [n+2,n+4])
               push_edge!(m, [n+4,n+3])
               push_edge!(m, [n+3,n+1])
           end
           m
end

#----for visualisation
function SemiAlgebraicTypes.mesh(rect::Rectangle)
           m = SemiAlgebraicTypes.mesh(Float64)
            push_vertex!(m, [rect.x_min, rect.y_min, 0.0])
            push_vertex!(m, [rect.x_max, rect.y_min, 0.0])
            push_vertex!(m, [rect.x_max, rect.y_max, 0.0])
            push_vertex!(m, [rect.x_min, rect.y_max, 0.0])
               push_edge!(m, [1,2])
               push_edge!(m, [2,3])
               push_edge!(m, [3,4])
               push_edge!(m, [4,1])
           m
end
