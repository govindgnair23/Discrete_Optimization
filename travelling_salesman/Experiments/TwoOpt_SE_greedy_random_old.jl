input_file = ARGS[1]

# Use greedy to get initial solution and then use two opt with simulated annealing
# Make a random change 
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
distance(city1::Tuple{Float64,Float64},city2::Tuple{Float64,Float64}) = sqrt((city1[1]-city2[1])^2 + (city1[2]-city2[2])^2)
distance(cities::Tuple{Int64,Int64}) = sqrt((coordinates[cities[1]].x-coordinates[cities[2]].x)^2 + (coordinates[cities[1]].y-coordinates[cities[2]].y)^2)



"""
Function to find city nearest to given city an distance provided a list of coordinates
Arguments:
    city: starting city
    coordinates: coordinates of all cities
    excluded_cities: cities that have been visited and hence should be excluded
"""
function get_nearest_city(city,coordinates,excluded_cities=[])
    #index i of this array gives distance from city to city i
    distances = Array{Float64}(undef,length(coordinates))
    #Set these to a high value so there are not selected
    distances[city] = Inf
    if !isempty(excluded_cities)
        distances[excluded_cities] .= Inf
    end

    for i ∈ range(1,stop=length(coordinates))
        if (city == i || i ∈ excluded_cities)
            continue
        end
        distances[i] = distance(coordinates[city],coordinates[i])
    end

    dist, closest_city = findmin(distances)

    return dist,closest_city
end


"""
Function to return a greedy optimal path given a set of co-ordinates
"""
function greedy_path(coordinates,seed=1)
    #Randomly select the starting city
    L = length(coordinates)
    starting_city = rand(1:L)
    shortest_path_ = Array{Int64}(undef,length(coordinates))
    shortest_path_[1] = starting_city
    total_distance = 0
    #Get nearest city
    for i ∈ 1:L-1
        current_city = shortest_path_[i]
        dist, closest_city = get_nearest_city(current_city,coordinates,shortest_path_[1:i])
        total_distance += dist
        shortest_path_[i+1] = closest_city
    end
    #Add distance to get back to starting city
    total_distance += distance(coordinates[shortest_path_[1]],coordinates[shortest_path_[end]])
    
    return total_distance,shortest_path_
end

"""
Function to run greedy path with random restarts N times and return the best path
"""
function run_greedy_N(N,coordinates)
    greedy_paths = [greedy_path(coordinates,i) for i in 1:N]
    best_path_i = argmin(greedy_paths)
    return greedy_paths[best_path_i]
end

#Function to eliminate out of bound error by getting element from requied index
function get_index(array,i)
    if i < 1
        return array[end]
    elseif i > length(array)
        return array[1]
    else
        return array[i]
    end
end

"""
Function to generate a new candidate with Two Opt algorithm and return new route distance
"""
function do_two_opt(route::Array{Int64,1},i::Int64,k::Int64,route_distance::Float64)

    #If i =1 and k = len(route) route is just reversed
    if i==1 && k==length(route)
        return route, route_distance 
    end
    L = length(route)
    new_route = Array{Int64,1}(undef,L)
    #Get nodes 1 to i-1
    if i > 1;new_route[1:i-1] = route[1:i-1];end
    
    #Take i to k and swap input_data
    new_route[i:k] = reverse(route[i:k])

    #Take k+1:n and add
    if k < L;new_route[k+1:end] = route[k+1:end];end

    replaced_edge_distance = distance(coordinates[get_index(route,i-1)],coordinates[route[i]]) + 
                                distance(coordinates[route[k]],coordinates[get_index(route,k+1)])
    new_edge_distance = distance(coordinates[get_index(route,i-1)],coordinates[new_route[i]]) + 
                        distance(coordinates[new_route[k]],coordinates[get_index(route,k+1)])

    Δ = new_edge_distance - replaced_edge_distance # Negative value means improvement
    return (route_distance + Δ,new_route)
end

"""
Function to calculate total distance of a route
"""
function calc_distance(route) 
    total_dist = 0

    for i in range(1,stop=length(route))
        total_dist += distance(coordinates[route[i]],coordinates[get_index(route,i+1)])
    end
    return total_dist
end

##Do exhaustive search using two opt
"""
Function to perform  local search using Metropolis Hastings
Inputs:
    route: current route
    route_distance: total distance of route
    T: temperature
"""
function local_search(route::Vector{Int64},route_distance::Float64,T::Float64)
    best_route = route
    best_route_distance = route_distance    
    N_cities = length(route)

    for i ∈ 1:N_cities-1
        for k ∈ i+1:N_cities
            new_route_distance,new_route = do_two_opt(route,i,k,route_distance)
            if new_route_distance < best_route_distance
                return new_route,new_route_distance
            else
                prob_metro = exp(-(new_route_distance - route_distance)/T)   #Calculate probability for metropolis evaluation
                if rand()< prob_metro
                    best_route = new_route
                    best_route_distance  = new_route_distance
                end
            end
        end
    end
    
    return best_route_distance,best_route
    
end


"""
Function to perform  local  random search using Metropolis Hastings
Inputs:
    route: current route
    route_distance: total distance of route
    T: temperature
"""
function local_search2(route::Vector{Int64},route_distance::Float64,T::Float64)
    best_route = route
    best_route_distance = route_distance    
    N_cities = length(route)

    i = rand(1:N_cities-1)
    k = rand(i+1:N_cities)

    new_route_distance,new_route, = do_two_opt(route,i,k,route_distance)
    if new_route_distance < best_route_distance
        return new_route_distance,new_route
    else
        prob_metro = exp(-(new_route_distance - route_distance)/T)   #Calculate probability for metropolis evaluation
        if rand()< prob_metro
            best_route = new_route
            best_route_distance  = new_route_distance
        end
    end
    
    return best_route_distance, best_route
    
end


## Do a random swap of two nodes to randomize solution
"""
Function that takes a route and swaps two nodes. Returns new route and route distance
Arguments
    route: route to be randomized
    route_distance: distance of route
    N_swaps: Number of random swaps to be carried out
"""
function random_swap(route,route_distance)
    L = length(route)
    node_indices = sample(1:L,2,replace=false) # nodes to be swapped
    discarded_edges = [((get_index(route,node_index-1),route[node_index]),
                        (route[node_index], get_index(route,node_index+1)))
                        for node_index ∈ node_indices]

    discarded_edge_distance = sum([distance(edge[1])+distance(edge[2]) for edge in discarded_edges])


    node1 = route[node_indices[1]]
    node2 = route[node_indices[2]]
    new_route = copy(route)
    new_route[node_indices[2]] = node1
    new_route[node_indices[1]] = node2

    new_edges = [((get_index(new_route,node_index-1),new_route[node_index]),
                        (new_route[node_index], get_index(new_route,node_index+1)))
                        for node_index ∈ node_indices]
    new_edge_distance =  sum([distance(edge[1])+distance(edge[2]) for edge in new_edges])

    return route_distance-discarded_edge_distance+new_edge_distance,new_route
    # new_route_distance = round(route_distance-discarded_edge_distance+new_edge_distance,digits=2)
    # new_route_distance_calc = round(calc_distance(new_route),digits=2)
    # return new_route_distance==new_route_distance_calc,new_route_distance,new_route_distance_calc
end

"""
Function to do multiple swaps for greater randomization
"""
function random_swap_N(route,route_distance,N::Int64)
    map(1:N) do _ 
         route_distance,route = random_swap(route,route_distance)
    end
    
    return route_distance,route
end


"""
Function to  do simulated annealing
T:: starting termperature
N_random:: Number of randomization to carry out
"""
function simulated_annealing(coordinates;budget::Int64,T::Float64,schedule::Float64,reset_point::Int64,
    N_random::Int64)
    best_route_distance, best_route = run_greedy_N(1,coordinates)
    temp = T
    counter = 0 # counter to keep track of number of cycles without improvement
    for i ∈ 1:budget
        #new_route,new_route_distance = local_search(best_route,best_route_distance,temp)
        new_route,new_route_distance = local_search(best_route,best_route_distance,temp)
        if new_route_distance < best_route_distance
            best_route = new_route
            best_route_distance = new_route_distance
            temp = schedule*temp
            counter = 0
        else
            counter +=1
        end

        if counter==reset_point
            best_route_distance, best_route =  random_swap_N(best_route,best_route_distance,N_random)
        end
    end

    return  best_route_distance,best_route 
end

##Do exhaustive search using two opt
"""
Function to perform exhaustive two opt, get best route and distance.
"""
function two_opt_exhaustive(route,route_distance) # To do: pass distance as an argument so as not to recompute route distance each time
    best_route = route
    best_route_distance = route_distance    
    N_cities = length(route)

    for i ∈ 1:N_cities-1
        for k ∈ i+1:N_cities
            new_route,new_route_distance = do_two_opt(route,i,k,route_distance)
            if new_route_distance < best_route_distance
                best_route = new_route
                best_route_distance  = new_route_distance
            end
        end
    end
    
    if best_route_distance == route_distance
        return  best_route_distance,best_route
    end
    
    return two_opt_exhaustive(best_route,best_route_distance)
end


#Use greedy to get initial solution
#initial_distance, initial_route = run_greedy_N(1,coordinates)

#optimal_route , optimal_distance =  two_opt_exhaustive(initial_route,initial_distance)
optimal_route,optimal_distance = simulated_annealing(coordinates;budget=10000,T=2.0,schedule=0.999,reset_point=100,N_random=10) #434
#Do an exhaustive two opt on result from simulated annealing
two_opt_exhaustive(optimal_route,optimal_distance)

#optimal_route , optimal_distance = simulated_annealing(coordinates,10000,0.8,0.99)
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
