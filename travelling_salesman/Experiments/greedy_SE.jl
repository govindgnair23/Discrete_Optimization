###### Greedy with simulated annealing########
### Solution here can be input to local search ####
input_file = ARGS[1]

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

using Random

#Function to calculate distance
distance(city1,city2) = sqrt((city1[1]-city2[1])^2 + (city1[2]-city2[2])^2)


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
    starting_city = rand(1:length(coordinates))
    shortest_path_ = []
    append!(shortest_path_,starting_city)
    total_distance = 0
    #Get nearest city
     while length(shortest_path_) < N_cities
        current_city = shortest_path_[end]
        dist, closest_city = get_nearest_city(current_city,coordinates,shortest_path_)
        total_distance += dist
        append!(shortest_path_,closest_city)
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

total_distance,solution = run_greedy_N(10,coordinates)
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