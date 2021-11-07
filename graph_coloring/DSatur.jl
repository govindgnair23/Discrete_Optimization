# Implementing the Degree of Saturation algorithm

input_file = ARGS[1]

f = open(input_file)
input_data = read(f,String)
close(f)

lines = split(input_data,'\n')

firstLine = split(lines[1])
global node_count = parse(Int64,firstLine[1])
global edge_count = parse(Int64,firstLine[2])

global edges = []
global nodes = []

for i ∈ range(2,stop=edge_count+1)
    line = lines[i]
    parts = split(line)
    edge = (parse(Int,parts[1]),parse(Int,parts[2]))
    append!(nodes,edge)
    push!(edges,edge)
end

using OffsetArrays
using StatsBase


"""
Function to get nodes that are connected to the given node
"""
function get_neighbors(edges,node)
    neighbors = []
    for edge ∈ edges
        if node ∈ edge
           push!(neighbors,setdiff(edge,node)[1])
        end
    end
    return neighbors
end


"""
Function to node rank by degrees of neighbors in descending order. 
Node with highest degree is the first element.
"""
function rank_nodes(nodes)
    degrees = countmap(nodes)
    degrees_ = sort(collect(degrees), by = x->x[2],rev=true)
    return [node_degree[1] for node_degree ∈ degrees_]
end


"""
Get most saturated nodes in a graph. In case of tie pick node with highest degree. If that is tied, randomize. <br>
Args:
    edges:: edges of the graph; global variable
    nodes:: Nodes in the graoh; global variable
    nodes_to_color:: Dictionary mapping node to color, assuming key is present only for colored nodes
    degrees:: Number of neighbors each node has
    ranked_nodes:: Nodes ranked by degrees to be used in case of ties
"""
function get_most_saturated(edges,nodes,ranked_nodes,nodes_to_color::Dict)
    nodes_ = unique(nodes)
    nodes_saturation = OffsetVector(ones(size(nodes_)).*-1,0:(node_count-1)) # Get saturation using node as index
    colored_nodes = collect(keys(nodes_to_color))
    delta_nodes = setdiff(nodes_,colored_nodes)
    #Get nodes that have been already colored
    for node ∈ delta_nodes
        neighbors = get_neighbors(edges,node)
        #Get  colors assigned to neighbors
        saturation = map(node->get(nodes_to_color,node,nothing),neighbors) # Get colors of neighboring nodes
        #Remove nothing
        saturation = filter(x-> !isnothing(x),saturation)
        #Get number of distinct colors and Assign that as saturation of node
        nodes_saturation[node] = length(unique(saturation))
    end
    #get most saturated node(s)
    most_saturated = findall(x-> x == maximum(nodes_saturation),nodes_saturation)
    #most_saturated = nodes[most_saturated_i.+1]
    if length(most_saturated) == 1
        return most_saturated[1]
    end
    
    # If there are ties get node with the highest degree
    # subgraph = get_subgraph(edges,nodes_to_color)
    # if !isempty(subgraph)
    #     nodes_ =   collect(Iterators.flatten(subgraph)) #Convert to list of nodes expected by rank_nodes function
    #     ranked_nodes_ = rank_nodes(nodes_)
    #     return ranked_nodes[1]
    # end

    #Subset just ties from list of ranked nodes
    ties = filter(x-> x ∈ most_saturated,ranked_nodes)
    #Return node with highest degree if there is a tie
    return ties[1]

end


"""
######## Redudnant function ##########
Function to return a subgraph(defined by edges only) induced by non-colored nodes.
Args:
    edges:: edges of the graph; global variable
    nodes_to_color:: Dictionary mapping node to color, assuming key is present only for colored nodes
Returns:
    edges without colored nodes
"""
function get_subgraph(edges,nodes_to_color)
    colored_nodes = collect(keys(nodes_to_color))
    #Retain edge only if it has non-colored nodes
    edges_ = [edge for edge in edges if (edge[1] ∉ colored_nodes && edge[2] ∉ colored_nodes)]
    return edges_
end

"""
Function to get number of colors required to color a graph and mapping of nodes to colors
"""
function graph_color_N(edges,nodes)
    N = length(unique(nodes)) #Get number of nodes in the graph
    color = 0
    nodes_to_color = Dict{Int64,Int64}()
    #Arrange vertices by decreasing order of degrees
    ranked_nodes = rank_nodes(nodes)
    #Choose vertex of maximal degree and assign color 0
    node = first(ranked_nodes)
    nodes_to_color[node] = color

    while length(nodes_to_color) < N
        #Choose a vertex with maximal saturation degree. If there is an equality, choose any vertex of maximal
        #degree 
        node_ = get_most_saturated(edges,nodes,ranked_nodes,nodes_to_color)
        #### Assign lowest possible color ###
        #Get neighbors
        neighbors_ = get_neighbors(edges,node_)
        #Get colors assigned to neighbors
        colors = [get(nodes_to_color,node_,nothing) for node_ in neighbors_]
        #Remove nothings
        colors = filter(x-> !isnothing(x),colors)

        max_color = maximum(colors)
        assigned_color = max_color+1
        

        #Choose lowest possible color  
        for col ∈ range(0,stop=max_color)
            if col ∉ colors
                assigned_color = col #Overwrite previously assigned color
            end
        end

        nodes_to_color[node_] = assigned_color
    end

    #Get number of colors used
    N_colors = maximum(values(nodes_to_color))+1

    return N_colors,nodes_to_color
end

N_colors,nodes_to_color = graph_color_N(edges,nodes)


"""
Function to render output in correct form
N_colors:: N_colors used for coloring the graph
nodes_to_color:: Dictionary mapping from node to color
node_count: Number of nodes in the graph
optimality_flag = 0
"""
function render_output(N_colors,nodes_to_color,node_count,optimality_flag)
    assignment = fill(0,node_count)
    for i in range(0,stop=node_count-1)
        assignment[i+1] = nodes_to_color[i]

    end

    println("$N_colors $optimality_flag")
    println(join(assignment," "))

end

render_output(N_colors,nodes_to_color,node_count,0)