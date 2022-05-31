##This is a direct replica of solution given in this tutorial 

using JuMP
import GLPK
import HiGHS
import Random
import Plots

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

dist_matrix = pairwise(distance,coordinates)
n_cities = length(coordinates)
"""
function to build tsp model
    d:: distance matrix
    n:: number of cities
"""
function build_tsp_model(d,n)
    model = Model(HiGHS.Optimizer)
    @variable(model, x[1:n,1:n],Bin,Symmetric)
    @objective(model,Min,sum(d .* x) /2)
    @constraint(model,[i in 1:n],sum(x[i,:]) == 2)
    @constraint(model,[i in 1:n], x[i,i] == 0)
    return model
end

"""
function to identify shortest subtour in current solution
"""
function subtour(edges::Vector{Tuple{Int,Int}},n)
    shortest_subtour,unvisited = collect(1:n),Set(collect(1:n))
    while !isempty(unvisited)
        this_cycle,neighbors = Int[],unvisited
        while !isempty(neighbors)
            current = pop!(neighbors)
            push!(this_cycle,current)
            if length(this_cycle) > 1
                pop!(unvisited, current)
            end
            neighbors = [j for (i,j) ∈ edges if i == current && j ∈ unvisited]
        end
        if length(this_cycle) < length(shortest_subtour)
            shortest_subtour = this_cycle
        end
    end
    return shortest_subtour
end


function selected_edges(x::Matrix{Float64},n)
    return Tuple{Int,Int}[(i,j) for i in 1:n, j in 1:n if x[i,j] > 0.5]
end

subtour(x::Matrix{Float64}) = subtour(selected_edges(x,size(x,1)),size(x,1))

subtour(x::AbstractMatrix{VariableRef}) = subtour(value.(x))


iterative_model = build_tsp_model(dist_matrix,n_cities)
set_optimizer_attribute(iterative_model,"log_to_console",false)
optimize!(iterative_model)
cycle = subtour(iterative_model[:x])
while 1 < length(cycle) < n_cities
    #println("Found cycle of length $(length(cycle))")
    S = [(i,j) for (i,j) in Iterators.product(cycle,cycle) if i < j]
    @constraint(
        iterative_model,
        sum(iterative_model[:x][i,j] for (i,j) in S) <= length(cycle)- 1,
    )
    optimize!(iterative_model)
    global cycle = subtour(iterative_model[:x])
end

optimal_distance = objective_value(iterative_model)


###################################
soln_mat = value.(iterative_model[:x])


#######Function to convert matrix solution to a route


function mat_to_route(soln_mat::Matrix{Float64},n_cities::Int64)
    prev_city = 1
    current_city = argmax(soln_mat[prev_city,:]) 
    optimal_route = Int64[prev_city,current_city]
    while length(optimal_route) < n_cities
        connected_cities = findall(x-> x==1.0,soln_mat[current_city,:])
        next_city = setdiff(connected_cities,optimal_route)[1]
        push!(optimal_route,next_city)
        current_city = next_city
    end
    return optimal_route
end

optimal_route = mat_to_route(soln_mat,n_cities)
optimal_route = optimal_route .- 1
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


