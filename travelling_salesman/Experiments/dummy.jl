## Dummy solution for problem #6 which seems intractable
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

#Function to calculate distance
distance(city1,city2) = sqrt((city1[1]-city2[1])^2 + (city1[2]-city2[2])^2)

solution = collect(0:N_cities-1)
#Function to conpute distance


function calc_distance()
    total_dist = 0

    for i in range(1,stop=N_cities-1)
        total_dist += distance(coordinates[i],coordinates[i+1])
    end
    total_dist += distance(coordinates[1],coordinates[end])
    return total_dist
end

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

total_dist = calc_distance()
render_output(total_dist,solution,0)