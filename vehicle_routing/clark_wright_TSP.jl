
using LinearAlgebra
using Plots
import HiGHS
using JuMP
using Distances
using Combinatorics
using Statistics
using StatsBase
#############Approach###############
#1) Use CLark Wright to assign cities to available vehicles
#2) Calculate centroid of cities assigned to each vehicle 
#3) Assign cities to these centroids - facility location problem
### 3.a) Try assigning unassigned city to closest centroid with excess capacity 
#4) Use TSP on each assigned city

######To do #####
##Objective function value is wrong. Error observed for some test cases.

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
N_customers = parse(Int64,firstLine[1]) #includes depot
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


function Base.:+(a::NamedTuple{(:x, :y),Tuple{Float64,Float64}}, b::NamedTuple{(:x, :y),Tuple{Float64,Float64}})
    return (x=a.x + b.x,y=a.y + b.y)
end

function Base.:/(a::NamedTuple{(:x, :y),Tuple{Float64,Float64}}, b::Int64)
    return (x=a.x/b,y=a.y/b)
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

#Function to calculate distance given coordinates
distance(city1,city2) = sqrt((city1.x-city2.x)^2 + (city1.y-city2.y)^2)

"""
Function to get distance matrix for a given array of cities + the depot
cites: cities to be included in the distance matrix
"""
function dist_matrix(coordinates,cities)
    cities_all = pushfirst!(coordinates[cities.+ 1],coordinates[1]) # Get depot
    return pairwise(distance,cities_all)
end

depot_coordinates = coordinates[1]
#Shift coordinates so depot is at origin
shifted_coordinates = [ coordinate - depot_coordinates for coordinate ∈ coordinates]
city_coordinates = shifted_coordinates[2:end] #cities to serve, excluding depot
#Function to calculate savings given a pair of coordinates
savings(city1,city2;depot=(x=0.0,y=0.0)) = distance(depot,city1)+distance(depot,city2) - distance(city1,city2)

###Get savings list###
combs = combinations(1:N_customers-1,2)
savings_dict = Dict()

for comb ∈ combs
    city1,city2 = city_coordinates[comb]
    savings_dict[comb] = savings(city1,city2)
end

savings_dict = sort(collect(savings_dict),by=x->x[2],rev=true)
#savings_matrix = pairwise(savings,shifted_coordinates[2:end])


"""
Check if a node is included in any of the routes and return route in which it is included
#is_in_routes(node,routes) = reduce(|,[node ∈ route for route in routes])
"""
function is_in_routes(node,routes)
    check = [node ∈ route for route in routes]
    return reduce(|,check), findfirst(check)
end

#is_not_interior(node,route) = node ∈ collect(route[2],route[end])
is_first(node,route) = node == route[2]#is node at non-interior point at begginnig of route
is_last(node,route) = node == route[end]#is node at non-interior point at end of route



"""
Function to merge two routes on two given cities
"""
function merge(route1,route2,city1,city2)
    route = nothing
    if is_first(city1,route1) && is_first(city2,route2) #If both cities come first in the route
        route = vcat([0],reverse(route2[2:end]),route1[2:end])
    elseif  is_first(city1,route1) && is_last(city2,route2)  #If one city comes first and second city comes last
        route = vcat([0],route2[2:end],route1[2:end])
    elseif is_last(city1,route1) && is_first(city2,route2)
        route = vcat([0],route1[2:end],route2[2:end])
    elseif  is_last(city1,route1) && is_last(city2,route2)
        route = vcat([0],route1[2:end],reverse(route2[2:end]))
    end
    return(route)
end


"""
Function to check if a city or multiple cities can be added to an existing route
"""
function check_capacity(route,cities;demand = demand[2:end] ,vehicle_capacity = vehicle_capacity)
    return sum(demand[route[2:end]]) + sum(demand[cities]) < vehicle_capacity ?  true  : false
end




"""
Given savings list as in clark wright algorithm, assign cities to vehicles
"""
function assign_to_vehicles(savings_dict;N_vehicles = N_vehicles)
    routes = [[0] for _ ∈ 1:N_vehicles]
    for (k,v) ∈ savings_dict
        next = findfirst(x->x==[0],routes) # Current empty route
        city1,city2 = k
        city1_in_route,route_no1 = is_in_routes(city1,routes)
        city2_in_route,route_no2 = is_in_routes(city2,routes)
    
     
        if !(city1_in_route|city2_in_route) # neither i nor j have already been assigned to a route, in which case a new route is initiated including both i and j. 
            !isnothing(next) && append!(routes[next],[city1,city2])
            #continue;
        end
    
    
    
        if city1_in_route && city2_in_route # if both are already part of some route
            route1 = routes[route_no1]
            route2 = routes[route_no2]
            if route_no1==route_no2;
                continue;#if both are in same route continue
            elseif check_capacity(route1,route2[2:end]) #merge routes if routes can be merged
                merged =  merge(route1,route2,city1,city2)
                if !isnothing(merged)
                    routes[route_no1] = merge(route1,route2 ,city1,city2)
                    routes[route_no2] = [0]
                end
            end
            continue;  
        end;
     
    
        #exactly one of the two points (i or j) has already been included in an existing route and that point is not interior to that route
        #(a point is interior to a route if it is not adjacent to the depot D in the order of traversal of points), in which case the link (i, j) is added to that same route. 
    
        if city1_in_route  # Then city 2 cannot be in the route
            route = routes[route_no1]
            if check_capacity(route,city2) # if there is room to add city to existing route
                if is_first(city1,route)
                    insert!(route,3,city2)
                elseif is_last(city1,route)
                    append!(route,city2)
                end
            elseif !isnothing(next) #there is a vacant route
                append!(routes[next],[city2])
            end
            continue;
        end
    
        if city2_in_route # Then city 1 cannot be in the route
            route = routes[route_no2]
            if check_capacity(route,city1) # if there is room to add city to existing route
                if is_first(city2,route)
                    insert!(route,3,city1)
                elseif is_last(city2,route)
                    append!(route,city1)
                end
            elseif !isnothing(next) #there is a vacant route
                append!(routes[next],[city1])
            end
            continue;
        end
    
        # If not assigned to any existing route, add to unassigned
    
    end

    return routes
end

routes  = assign_to_vehicles(savings_dict);

"""
Function to convert routes into an array of assignments
"""
function routes_to_assignment(routes;N_customers=N_customers)
    assignment = zeros(N_customers-1)
    for (i,route) ∈ enumerate(routes)
        assignment[route[2:end]] .= i
    end
    return assignment
end


assignment = routes_to_assignment(routes)

######Determine if there are any unassigned cities, if yes identify cluster centroids and reassign cities ####
all_assigned_flag = isempty(findall(x->x==0,assignment)) #True if all cities have been assigned to routes 
 

"""
Function to get centroid for each route using assigned cities##
"""
function centroids(routes::Array{Array{Int64,1}},coordinates::Array{Any,1})
    centroids = []
    for route ∈ routes
        push!(centroids,sum(coordinates[route.+1])/length(route))
    end

    return centroids
end



"""
Function to build facility location model
"""
function build_fl_model(n_facilites,n_cities,a,q,c,f)
    #Model
    model = Model(HiGHS.Optimizer)

    #Variables
    @variable(model,y[1:n_facilites],Bin);
    @variable(model, x[1:n_cities,1:n_facilites], Bin);

    #Each city is assigned only to 1 facility
    @constraint(model,customer_service[i in 1:n_cities],sum(x[i,:])==1);
    #Capacity constraint
    @constraint(model, capacity,x'a .<= (q.*y))

    #objective
    @objective(model, Min, f'y + sum( c .* x));

    return model
    
end

"""
Function to assign cities to facilties located in centroids
"""
function assign_to_facilities(routes,coordinates;N_vehicles=N_vehicles,N_customers=N_customers,demand=demand)
    facilities = centroids(routes,coordinates)
    n_facilites = N_vehicles    
    
    cities = coordinates[2:end]
    n_cities = N_customers - 1
    @assert length(cities) == n_cities

    ###### Set up facility location problem #########
    #Fixed costs for a facility: assumed to be 0#
    f = zeros(n_facilites)
    
    #Demands from each city
    a = convert(Array{Int64},demand[2:end])

    
    #Capacity of depot/vehicle
    q = repeat([vehicle_capacity],N_vehicles)

    
    #Transportation costs/Cost of assigning a city to a route
    c = zeros(n_cities,n_facilites)
    for i ∈ 1:n_cities
        for j ∈ 1:n_facilites 
            c[i,j] = distance(cities[i],facilities[j])
        end
    end

    #Model
    cfl = build_fl_model(n_facilites,n_cities,a,q,c,f)
    set_optimizer_attribute(cfl,"log_to_console",false)
    optimize!(cfl)

    
    x_ = round.(value.(cfl[:x]),digits=1)

    assignment = map(eachrow(x_)) do y #Error here
        findall(x->x >=1.0,y)[1]
    end

    return assignment
end

if !(all_assigned_flag)
    assignment = assign_to_facilities(routes,coordinates)
end


###############Set up MIP TSP to find optimal route for each vehicle###


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
function vehicle_routing(vehicle_assignment,coordinates;N_vehicles=N_vehicles)
    total_distance = 0
    routes = [] 
    for vehicle ∈ 1:N_vehicles 
        cities = findall(x->x==vehicle,vehicle_assignment)
        ##If only one city assigned to a vehicle, no need to do TSP
        if length(cities) == 1
            push!(routes,[0,cities[1]])
            continue
        end

        # if no cities are assigned to the vehicle, treat appropriately
        if isempty(cities)
            push!(routes,[0])
            continue
        end



        d_matrix = dist_matrix(coordinates,cities)
        n_cities = length(cities)+1 

        iterative_model  = build_tsp_model(d_matrix,n_cities)
        set_optimizer_attribute(iterative_model,"log_to_console",false)
        optimize!(iterative_model)
        cycle = subtour(iterative_model[:x])

        while 1 < length(cycle) < n_cities
            S = [(i,j) for (i,j) in Iterators.product(cycle,cycle) if i < j]
            @constraint(
                iterative_model,
                sum(iterative_model[:x][i, j] for (i, j) in S) <= length(cycle) - 1,
            )
            optimize!(iterative_model)
            cycle = subtour(iterative_model[:x])
        end

        soln_mat = value.(iterative_model[:x])
        optimal_route = mat_to_route(soln_mat,n_cities)
        pushfirst!(cities,0) # Add depot to route
        #Order cities as per route
        ordered_cities = cities[optimal_route]
        #Add distance to return to depot
        total_distance += objective_value(iterative_model) + distance(coordinates[1],coordinates[ordered_cities[end]])
        push!(routes,ordered_cities)
    end

    return routes,total_distance
end

routes,total_distance = vehicle_routing(assignment,coordinates)




#######Plot the routes obtained####
using Plots

"""
Function to plot routes
"""
function plot_routes(routes)
    p = plot()
    for route ∈ routes
        push!(route,0)
        x=[coord.x for coord in coordinates[route.+1]]
        y=[coord.y for coord in coordinates[route.+1]]
        plot!(x,y)
    end
    p
end


#plot_routes(routes)

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