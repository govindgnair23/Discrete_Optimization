input_file = ARGS[1]

f = open(input_file)
input_data = read(f,String)
close(f)

lines = split(input_data,'\n')

global N_cities = parse(Int64,lines[1])
coordinates = []

for i ∈ 2:N_cities+1
    line = lines[i]
    coords = split(line)
    push!(coordinates,(x=parse(Float64,coords[1]),y=parse(Float64,coords[2])))
end

using LinearAlgebra

#Declare a distance matrix
Dist = Array{Float64,2}(undef,N_cities,N_cities)

#Function to calculate distance
distance(city1,city2) = sqrt((city1[1]-city2[1])^2 + (city1[2]-city2[2])^2)

for i ∈ 1:N_cities
    for j ∈ 1:N_cities
        Dist[i,j] = distance(coordinates[i],coordinates[j])
    end
end

#Assign high value to diagonal elements
Dist[diagind(Dist)] .= Inf

"""
Function to find shortest path given a distance matrix
Arguments:
    D:: Distance Matrix
"""
function shortest_path(D::Array{Float64,2})
    current_city = 1  #Assume we start at city 1
    covered_cities = [1]
    total_distance = 0 #distance covered
    while length(covered_cities) < N_cities
        #Won't have to visit current city again so set distances to current city to a very high value
        D[:,current_city].= Inf
        #Get distances to citites that can be visited from current city
        distances = D[current_city,:]
        #Get closest city
        dist, closest_city = findmin(distances)
        append!(covered_cities,closest_city)
        total_distance += dist
        current_city = closest_city
    end

    return total_distance,covered_cities
end

total_distance,solution = shortest_path(Dist)
#Also need to add distance from last city to first city
total_distance += Dist[1,last(solution)]

#Solution expects first city to be 0, so deduct 1
solution = solution .- 1

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

render_output(total_distance,solution,0)