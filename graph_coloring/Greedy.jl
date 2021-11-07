#Implementing the Welsh Powell Graph Coloring Algorithms
input_file = ARGS[1]

f = open(input_file)
input_data = read(f,String)
close(f)

lines = split(input_data,'\n')

firstLine = split(lines[1])
node_count = parse(Int64,firstLine[1])
edge_count = parse(Int64,firstLine[2])

global edges = []
global nodes = []

for i ∈ range(2,stop=edge_count+1)
    line = lines[i]
    parts = split(line)
    edge = (parse(Int,parts[1]),parse(Int,parts[2]))
    append!(nodes,edge)
    push!(edges,edge)
end

using StatsBase



"""
Function to node rank by degrees of neighbors in descending order. 
Node with highest degree is the first element.
"""
function rank_nodes(nodes)
    degrees = countmap(nodes)
    degrees_ = sort(collect(degrees), by = x->x[2],rev=true)
    return [edge[1] for edge ∈ degrees_]
end


"""
Function to get nodes that are connected to the given node
"""
function get_neighbors(edges,node)
    neighbors = Set{Int64}()
    for edge ∈ edges
        if node ∈ edge
           push!(neighbors,setdiff(edge,node)[1])
        end
    end
    return neighbors
end


"""
Function to check if node2 is a neighbor of node1
"""
function is_neighbor(node1,node2)
    return (node1,node2) ∈ edges
end


"""
Function to get number of colors required to color a graph
"""
function graph_color_N(edges,nodes)
    ranked_nodes = rank_nodes(nodes)
    color_to_nodes = Dict()
    color = 0
    while !isempty(ranked_nodes)
        node = popfirst!(ranked_nodes) #Get higest ranked node
        nodes_to_color = [node]
        excluded_nodes = Set(get_neighbors(edges,node))
        for node_ ∈ ranked_nodes #loop through  remaining nodes
            if node_ ∉ excluded_nodes
                push!(nodes_to_color,node_)
                union!(excluded_nodes,get_neighbors(edges,node_))
            end
        end
        color_to_nodes[color] = nodes_to_color #Map nodes to be colored to a specific color..
        setdiff!(ranked_nodes,nodes_to_color) #All non neighbors have been assigned a color so remove
        color+=1
    end
    return color,color_to_nodes
end

N_colors,color_to_nodes = graph_color_N(edges,nodes)

"""
Function to render output in correct form
item_count: No of items to select from
selected_items: List of selected items by index
value: Value of solution
optimality_flag = 0
"""
function render_output(N_colors,color_to_nodes,node_count,optimality_flag)
    assignment = fill(0,node_count)
    for i in range(1,stop=N_colors-1)
        #Get nodes assigned with color i
        assigned_nodes = color_to_nodes[i]
        assignment[assigned_nodes.+1].= i

    end

    println("$N_colors $optimality_flag")
    println(join(assignment," "))

end

render_output(N_colors,color_to_nodes,node_count,0)
