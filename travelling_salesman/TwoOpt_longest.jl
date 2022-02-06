input_file = ARGS[1]

# Implement the 2-Opt technique but with a more limited neigborhood serach.
# Swap only the two longest edges
f = open(input_file)
input_data = read(f,String)
close(f)

lines = split(input_data,'\n')

global N_cities = parse(Int64,lines[1])
global coordinates = []

for i ∈ 2:N_cities+1
    line = lines[i]
    coords = split(line)
    push!(coordinates,(x=parse(Float64,coords[1]),y=parse(Float64,coords[2])))
end

using StatsBase
using Random

#Function to calculate distance given coordinates
distance(city1,city2) = sqrt((city1[1]-city2[1])^2 + (city1[2]-city2[2])^2)

#Function to eliminate out of bound error by getting element from requied index
function get_index(array::Vector{Int64},i::Int64)
    if i < 1
        return array[end]
    elseif i > length(array)
        return array[1]
    else
        return array[i]
    end
end

"""
Function to generate a new candidate with Two Opt algorithm and return new route distance.
Reverse segment of route from i+1 to k (inclusive)
"""
function do_two_opt(route::Array{Int64,1},i::Int64,k::Int64,route_distance::Float64,L::Int64)

    #If i =1 and k = len(route) route is just reversed
    # if i==1 && k==length(route)
    #     return route, route_distance 
    # end

    new_route = similar(route)
    #Get nodes 1 to i-1
    # if i > 1;new_route[1:i-1] = route[1:i-1];end
    
    #Retain 1:i
    new_route[1:i] = route[1:i]

    #Take i+1 to k and swap input_data
    new_route[i+1:k] = reverse(route[i+1:k])

    #Take k+1:n and add
    if k < L;new_route[k+1:end] = route[k+1:end];end

    replaced_edge_distance = distance(coordinates[route[i]],coordinates[get_index(route,i+1)]) + 
                                distance(coordinates[route[k]],coordinates[get_index(route,k+1)])
    new_edge_distance = distance(coordinates[new_route[i]],coordinates[get_index(new_route,i+1)]) + 
                                distance(coordinates[route[k]],coordinates[get_index(route,k+1)])

    Δ = new_edge_distance - replaced_edge_distance # Negative value means improvement
    return (new_route,route_distance + Δ)
end

#Function to compute distance traversed in a given route and return distance along with edges in
# descending order of length
# N:: number of (longest)edges to be returned 
##Reuturns length of route and starting index of nodes with longest edges
## If function return 2&4 , route[2]-route[3] and route[4]-route[5] are longest


"""
Function to calculate distance of a route along with N longest edges
Arguments:
    route:: route whose distance is to be calculated
    N:: Number of longest edges
Returns:
    Length of route
    Starting index of node with longest edges. If function return 2&4 , route[2]-route[3] and route[4]-route[5] are longest
"""
#hash_table = Dict() #Keep track of distances to avoid recomputing

function calc_distance(route,N) 
    total_dist = 0
    edge_lengths = []

    for i in range(1,stop=length(route))
        edge_length = distance(coordinates[route[i]],coordinates[get_index(route,i+1)])
        push!(edge_lengths,edge_length)
        total_dist += edge_length
    end
    longest_edges = sortperm(edge_lengths,rev=true)
    return total_dist, sort(longest_edges[1:N])
end


##Do limited two opt search by replacing longest edges
"""
Function to perform two opt, get best route and distance.
Arguements::
    route:: route from which to begin local search
    N: number of nodes to consider as candidates for two opt swap
"""
function two_opt_longest(route,N) 
    route_distance,swap_candidates = calc_distance(route,N)
    best_route = route
    best_route_distance = route_distance    
    N_cities = length(route)


    for i ∈ swap_candidates[1:end-1]
        for k ∈ swap_candidates[2:end]
            new_route,new_route_distance = do_two_opt(route,i,k,route_distance,N_cities)
            if new_route_distance < best_route_distance
                best_route = new_route
                best_route_distance  = new_route_distance
            end
        end
    end
    
    if best_route_distance == route_distance
        return  best_route_distance, best_route
    end
    
    return two_opt_longest(best_route,N)

    
end



if N_cities <10000
    optimal_distance,optimal_route =  two_opt_longest(shuffle!(collect(1:N_cities)),10)
else
    optimal_distance,optimal_route =  two_opt_longest(shuffle!(collect(1:N_cities)),2)
end
optimal_route = optimal_route .- 1

#############################################################################
"""
Function to render output in correct form
Arguments: 
    distance: Total distance traversed by salesman
    solution: sequence of cities to visit
    optimality_flag = 0
"""
function render_output(distance,solution,optimality_flag=0)
    println("$distance $optimality_flag")
    println(join(solution," "))
end

render_output(optimal_distance,optimal_route,0)

## To do: This works only if starting and ending node are same as in optimal route

# i and k can range from 2 to L -1 where K>i
#Edges replaced are i-1:i and k:k+1


#Need to update above. should be able to set i to 1
#A-B-C-D -> B-A-C-D i=1,k=2
#A-B-C-D -> A-B-D-C ,i=3,k=4


######################### Testing ##############

## TO DO: Create a Struct with wrap around indexing


# route = collect(1:5)
# calc_distance(route)

# new_route,delta = do_two_opt(route,2,4)
# new_route,delta = do_two_opt(route,2,3)

# calc_distance(new_route)


#optimal_distance,optimal_route =  two_opt_longest(shuffle!(collect(1:N_cities)),200)