
using LinearAlgebra
using Plots
import HiGHS
using JuMP
using Distances

input_file = ARGS[1]

file = open(input_file)
input_data = read(file,String)
close(file)

lines = split(input_data,'\n')

# demands
global demand = []
# coordinates
global coordinates = Any[]


firstLine = split(lines[1])
#Number of facilites
N_customers = parse(Int64,firstLine[1]) 
N_vehicles = parse(Int64,firstLine[2])
vehicle_capacity =  parse(Int64,firstLine[3])

for i ∈ range(2,stop=N_customers+1)
    line = lines[i]
    parts = split(line)
    push!(demand,parse(Float64,parts[1]))
    push!(coordinates,(x = parse(Float64,parts[2]), y = parse(Float64,parts[3])))
end


#Function  subtract named tuples
function Base.:-(a::NamedTuple{(:x, :y),Tuple{Float64,Float64}}, b::NamedTuple{(:x, :y),Tuple{Float64,Float64}})
    return (x=a.x - b.x,y=a.y -b.y)
end


########Convert to coordinates to polar coordinates
function to_polar((x,y))
   if x > 0 && y >= 0 #First quadrant
        theta = atan(y/x) 
   elseif x <= 0 && y > 0 #Second quadrant
        theta = π - atan(y/abs(x))
   elseif x < 0 && y <= 0 # Third quadrant
        theta = π + atan(y/x)
   else #Fourth quadrant
        theta = 2π - atan(abs(y)/x)
   end
   return (sqrt(x^2+y^2),theta)
end

"""
Assign customers to vehicles using the petal/sweep algorithm.
"""
function assign_vehicles(demand, coordinates,N_customers;chart=false)
    #Depot has no demand
    demand = demand[2:end]
    
    #First coordinate is the depot
    depot_coordinates = coordinates[1]
    #Shift coordinates so that depot is at the origin
    shifted_coordinates = [ coordinate - depot_coordinates for coordinate ∈ coordinates]
    polar_coordinates  = map(to_polar,shifted_coordinates[2:end])

    #Get theta values for each customer 
    theta = [x[2] for x in polar_coordinates]

    #Get sequence of customer in increasing order of theta values
    order_customers = sortperm(theta)

    #Sequence demand in above order
    demand_ = demand[order_customers]

    #Get cumulative demand
    cum_demand = cumsum(demand_)

    customer_counter = 0
    vehicle_counter = 1
    c = 1
    vehicle_assignment = Dict()
    while customer_counter < N_customers - 1 # depot 0 is not a customer
        c_old = c
        #Get rank of customer that causes vehicle to fill out
        c =  findfirst( x-> x > vehicle_counter*vehicle_capacity,cum_demand) 

        if isnothing(c)
            c = N_customers
        end
        #Customer from current value of customer_counter to c belong to vehicle number equal to vehicle vehicle_counter
        customer_counter += c - c_old 
        assigned_customers = order_customers[c_old:c-1]
        #Assign the customers to first vehicle
        vehicle_assignment[vehicle_counter] = assigned_customers
        vehicle_counter+=1
    end

    if chart == true
        ##Convert array of named tuples to simple arrays
        x,y = [getindex.(shifted_coordinates[2:end],i) for i in 1:2]
        plot()
        for k ∈ keys(vehicle_assignment)
            v = vehicle_assignment[k]
            plot!(x[v],y[v],color = k,seriestype = :scatter,legend = false)
        end
        display(current())
    end

    return vehicle_assignment

end

vehicle_assignment = assign_vehicles(demand,coordinates,N_customers)


###############Set up MIP TSP to find optimal route for each vehicle###
#Function to calculate distance given coordinates
distance(city1,city2) = sqrt((city1.x-city2.x)^2 + (city1.y-city2.y)^2)

"""
Function to get distance matrix for a given array of cities + the depot
"""
function dist_matrix(coordinates,cities)
    cities_all = push!(coordinates[cities.+ 1],coordinates[1]) # Get depot
    return pairwise(distance,cities_all)
end



"""
Function to build a TSP Model with MIP
"""

function build_tsp_model(d,n)
    model = Model(HiGHS.Optimizer)
    @variable(model,x[1:n,1:n],Bin,Symmetric)
    @objective(model,Min, sum(d .* x)/2)
    @constraint(model,[ i in 1:n],sum(x[i,:])==2)
    @constraint(model,[i in 1:n], x[i,i] == 0)
    return model
end


function subtour(edges::Vector{Tuple{Int,Int}},n)
    shortest_subtour , unvisited = collect(1:n), Set(collect(1:n))
    while !isempty(unvisited)
        this_cycle , neighbors = Int[], unvisited
        while !isempty(neighbors)
            current = pop!(neighbors)
            push!(this_cycle, current)
            if length(this_cycle) > 1
                pop!(unvisited,current)
            end

            neighbors = [j for (i,j) in edges if i == current && j in unvisited]
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


### Solve one sub-problem

d_matrix = dist_matrix(coordinates,vehicle_assignment[1])
n_cities = length(vehicle_assignment[1])+1 


###Solve iteratively
iterative_model  = build_tsp_model(d_matrix,n_cities)
optimize!(iterative_model)
cycle = subtour(iterative_model[:x])

while 1 < length(cycle) < n_cities
    S = [(i,j) for (i,j) in Iterators.product(cycle,cycle) if i < j]
    @constraint(
        iterative_model,
        sum(iterative_model[:x][i, j] for (i, j) in S) <= length(cycle) - 1,
    )
    optimize!(iterative_model)
end

soln_mat = value.(iterative_model[:x])

mat_to_route(soln_mat,5)
###################################
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

##Function to solve full problem###
function vehicle_routing(vehicle_assignment,coordinates)
    total_distance = 0
    routes = []
    vehicles =  sort(collect(keys(vehicle_assignment)))
    for vehicle ∈ vehicles
        d_matrix = dist_matrix(coordinates,vehicle_assignment[vehicle])
        n_cities = length(vehicle_assignment[vehicle])+1 

        iterative_model  = build_tsp_model(d_matrix,n_cities)
        optimize!(iterative_model)
        cycle = subtour(iterative_model[:x])

        while 1 < length(cycle) < n_cities
            S = [(i,j) for (i,j) in Iterators.product(cycle,cycle) if i < j]
            @constraint(
                iterative_model,
                sum(iterative_model[:x][i, j] for (i, j) in S) <= length(cycle) - 1,
            )
            optimize!(iterative_model)
        end

        soln_mat = value.(iterative_model[:x])
        route_order = mat_to_route(soln_mat::Matrix{Float64},n_cities::Int64)
        cities  = copy(vehicle_assignment[vehicle])
        pushfirst!(cities,0)
        optimal_route  = cities[route_order]
        total_distance += objective_value(iterative_model)
        push!(routes,optimal_route)
    end

    return routes,total_distance
end

routes,total_distance = vehicle_routing(vehicle_assignment,coordinates)

"""
Function to render output in correct form
Arguments: 
    distance: Total distance traversed by all vehicles
    solution:  list of lists where each sublist is route
    optimality_flag = 0
"""
function render_output(distance,solution,N_vehicles;optimality_flag=0)
    println("$distance $optimality_flag")
    N_assigned_vehicles = length(solution)
    for i ∈ 1:N_assigned_vehicles
        println(join(solution[i]," ")," 0")
    end
    if N_vehicles >N_assigned_vehicles
        for i ∈ N_assigned_vehicles+1:N_vehicles
            println("0 0")
        end
    end
end

render_output(total_distance,routes,N_vehicles)