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

#Function to calculate length of two edges in a sequence
#Computes distance node1-node2 and node3-node4
function distance2opt((node1,node2,node3,node4),)
     return distance(coordinates[node1],coordinates[node2]) + distance(coordinates[node3],coordinates[node4])
end

"""
Function to generate a new candidate with Two Opt algorithm and return Δ = dist(new_route) - dist(old_route)
"""
function do_two_opt(route::Array{Int64,1},i::Int64,k::Int64)
    L = length(route)
    new_route = Array{Int64,1}(undef,L)
    #Get nodes 1 to i-1
    #append!(new_route,(route[1:i-1]))
    new_route[1:i-1] = route[1:i-1]
    #Take i to k and swap input_data
    #append!(new_route,(reverse(route[i:k])))
    new_route[i:k] = reverse(route[i:k])
    #Take k+1:n and add
    #append!(new_route,(route[k+1:end]))
    new_route[k+1:end] = route[k+1:end]

    replaced_edge_distance = distance(coordinates[route[i-1]],coordinates[route[i]]) + 
                             distance(coordinates[route[k]],coordinates[route[k+1]])
    new_edge_distance = distance(coordinates[new_route[i-1]],coordinates[new_route[i]]) + 
                        distance(coordinates[new_route[k]],coordinates[new_route[k+1]])

    Δ = new_edge_distance - replaced_edge_distance # Negative value means improvement
    return (new_route,Δ)
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

##Do exhaustive search using two opt
"""
Function to perform two opt, get best route and distance.
"""
function two_opt_exhaustive(route) #Current implementation requires to cycle through all ordering representations
                                   # E.g. 4,3,1,5,2 and 3.1.5.2.4
    best_route = route
    for i ∈ 2:N_cities-2
        for k ∈ i+1:N_cities-1
            new_route,Δ = do_two_opt(best_route,i,k)
            if Δ < 0
                return two_opt_exhaustive(new_route)
            end  
        end
    end
    best_route_distance = calc_distance(best_route)
    return best_route_distance, best_route

end


two_opt_exhaustive(collect([3,4,2,1,5]))

## To do: This works only if starting and ending node are same as in optimal route

# i and k can range from 2 to L -1 where K>i
#Edges replaced are i-1:i and k:k+1


#Need to update above. should be able to set i to 1
#A-B-C-D -> B-A-C-D i=1,k=2
#A-B-C-D -> A-B-D-C ,i=3,k=4


######################### Testing ##############




route = collect(1:5)
calc_distance(route)

new_route,delta = do_two_opt(route,2,4)
new_route,delta = do_two_opt(route,2,3)

calc_distance(new_route)