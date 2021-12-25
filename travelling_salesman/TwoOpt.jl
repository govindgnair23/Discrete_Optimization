# Implement the 2-Opt technique with simulated annealing
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

#Function to calculate distance given coordinates
distance(city1,city2) = sqrt((city1[1]-city2[1])^2 + (city1[2]-city2[2])^2)

#Function to eliminate out of bound error by getting element from requied index
function get_index(array,i)
    if i < 1
        return array[end]
    elseif i > N_cities
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
    return (new_route,route_distance + Δ)
end

#Function to compute distance traversed in a given route..
function calc_distance(route) # Not working correctly.
    total_dist = 0
    for i in range(1,stop=length(route)-1)
        total_dist += distance(coordinates[route[i]],coordinates[route[i+1]])
    end
    total_dist += distance(coordinates[route[1]],coordinates[route[end]])
    return total_dist
end

function check_improvement!(i,k,best_route,best_route_distance,improved)
    if i == 1 & k == N_cities
        return best_route, best_route_distance, improved
    end
    new_route,Δ = do_two_opt(best_route,i,k)
    if Δ < 0
        best_route = new_route
        best_route_distance += Δ
        improved = true
    end
    #return best_route, best_route_distance, improved
end



##Do exhaustive search using two opt
"""
Function to perform two opt, get best route and distance.
"""
function two_opt_exhaustive(route) 
    route_distance = calc_distance(route)
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
        return  best_route_distance, best_route
    end
    
    return two_opt_exhaustive(best_route)

    
end




two_opt_exhaustive(collect([3,2,4,1,5]))

## To do: This works only if starting and ending node are same as in optimal route

# i and k can range from 2 to L -1 where K>i
#Edges replaced are i-1:i and k:k+1


#Need to update above. should be able to set i to 1
#A-B-C-D -> B-A-C-D i=1,k=2
#A-B-C-D -> A-B-D-C ,i=3,k=4


######################### Testing ##############

## TO DO: Create a Struct with wrap around indexing


route = collect(1:5)
calc_distance(route)

new_route,delta = do_two_opt(route,2,4)
new_route,delta = do_two_opt(route,2,3)

calc_distance(new_route)

