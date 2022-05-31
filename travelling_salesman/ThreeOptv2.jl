input_file = ARGS[1]

# Implement exhaustive three opt technique
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
Function to reverse the segment inclusive of given indices
"""
function reverse_segment(route::Array{Int64,1},i::Int64,j::Int64,coordinates)
    L = length(route)
    new_route = copy(route)

    if j > L
        i = 1
        j = L
    end 
    new_route[i:j] = reverse(route[i:j]) 
    
    replaced_edge_distance = distance(coordinates[get_index(route,i-1)],coordinates[route[i]]) + 
                                distance(coordinates[route[j]],coordinates[get_index(route,j+1)])
    new_edge_distance = distance(coordinates[get_index(new_route,i-1)],coordinates[new_route[i]]) + 
                        distance(coordinates[new_route[j]],coordinates[get_index(new_route,j+1)])

    Δ = new_edge_distance - replaced_edge_distance #negative means reduction in distance

    return new_route,Δ
end


"""
Function to apply three opt move.
Returns: new route and distance of new route

"""
function do_three_opt_move(route::Array{Int64,1},i::Int64,j::Int64,k::Int64,coordinates,route_distance)
   
    L = length(route)
    case1,Δ1 = reverse_segment(route,i+1,j,coordinates)#a'bc
    case2,Δ2 = reverse_segment(route,k+1,L,coordinates)#abc'
    case3,Δ3 = reverse_segment(route,j+1,k,coordinates)#ab'c
    case4,Δ4 = reverse_segment(case3,k+1,L,coordinates)#ab'c'
    Δ4 += Δ3
    case5,Δ5 = reverse_segment(case3,i+1,j,coordinates)#a'b'c
    Δ5 += Δ3
    case6,Δ6 = reverse_segment(case1,k+1,L,coordinates)#a'bc'
    Δ6 += Δ1
    case7,Δ7 = reverse_segment(case4,i+1,j,coordinates)#a'b'c'
    Δ7 += Δ4

    Cases = [case1,case2,case3,case4,case5,case6,case7]
    Delta = [Δ1,Δ2,Δ3,Δ4,Δ5,Δ6,Δ7]

    min_ = argmin(Delta)
    if Delta[min_] < 0
        return Cases[min_],route_distance+ Delta[min_]
    end

    return route, route_distance


end

"""
Function to calculate length of route and return N longest edges
"""
function calc_distance(route,N) 
    total_dist = 0
    edge_lengths = []

    for i in range(1,stop=length(route))
        edge_length = distance(coordinates[route[i]],coordinates[get_index(route,i+1)])
        push!(edge_lengths,edge_length)
        total_dist += edge_length
    end
    longest_edges = sortperm(edge_lengths,rev=true)
    return total_dist, longest_edges[1:N] #Gives indices of starting node
                                          #E.g. if longest_edge  = 3 => route[3] - route[4] is the longest edge
end

"""
Function to calculate length of route  only
"""
function calc_distance(route) 
    total_dist = 0

    for i in range(1,stop=length(route))
        edge_length = distance(coordinates[route[i]],coordinates[get_index(route,i+1)])
        total_dist += edge_length
    end

    return total_dist
end


##Do limited three opt search by replacing longest edges
"""
Function to perform two opt, get best route and distance.
Arguements::
    route:: route from which to begin local search
    N: number of edges to consider for three opt move
"""
function three_opt_longest(route,N,coordinates) 
    route_distance,edge_candidates = calc_distance(route,N)
    edge_candidates = setdiff(edge_candidates,length(route))
    best_route = route
    best_route_distance = route_distance    
    N_ =  length(edge_candidates)
    for i ∈ range(1,stop=N_-2)
        i_,j_,k_ =  sort(edge_candidates[i:i+2])
        new_route,new_route_distance = do_three_opt_move(route,i_,j_,k_,coordinates,route_distance)
        if new_route_distance < best_route_distance
            best_route = new_route
            best_route_distance  = new_route_distance
        end
    end

    if best_route_distance == route_distance
        return  best_route_distance, best_route
    end
    
    return three_opt_longest(best_route,N,coordinates) 


end

"""
Function to do three opt on random edges. N specifies number of times to carry out this operation
"""
function three_opt_random(route,route_distance,N,coordinates;N_cities=N_cities) 
    best_route = route
    best_route_distance = route_distance    
    for i ∈ range(1,stop=N)
        i_,j_,k_ =  sort(sample(1:N_cities-1,3,replace=false))
        new_route,new_route_distance = do_three_opt_move(best_route,i_,j_,k_,coordinates,best_route_distance)
        if new_route_distance < best_route_distance
            best_route = new_route
            best_route_distance  = new_route_distance
        end
    end

    return  best_route_distance, best_route

end

"""
Function to do three opt on random edges with simulated_annealing. 
"""
function three_opt_se(route,route_distance,coordinates,T;N_cities=N_cities) 
    best_route = route
    best_route_distance = route_distance    
   
    i_,j_,k_ =  sort(sample(1:N_cities-1,3,replace=false))
    new_route,new_route_distance = do_three_opt_move(best_route,i_,j_,k_,coordinates,best_route_distance)
    if new_route_distance < best_route_distance
        best_route = new_route
        best_route_distance  = new_route_distance
    else
        prob_metro = exp(-(new_route_distance - best_route_distance)/T)   #Calculate probability for metropolis evaluation
        if rand()< prob_metro
            best_route = new_route
            best_route_distance  = new_route_distance
        end
    end
 

    return  best_route_distance, best_route

end


"""
Function to  do simulated annealing
T:: starting termperature
N_random:: Number of randomization to carry out
reset_point:: Number of cycles with no improvement before randomizing solution
"""
function simulated_annealing(coordinates;budget::Int64,T::Float64,schedule::Float64,reset_point::Int64,N_cities=N_cities)
    best_route = shuffle(collect(1:N_cities))
    best_route_distance = calc_distance(best_route)
    intense_route_distance , intense_route = copy(best_route_distance),copy(best_route)
    temp = T
    counter = 0 # counter to keep track of number of cycles without improvement
    for i ∈ 1:budget
        new_route_distance,new_route = three_opt_se(best_route,best_route_distance,coordinates,T)
        if new_route_distance < best_route_distance
            best_route = new_route
            best_route_distance = new_route_distance
            temp = schedule*temp
            counter = 0 # if new optimum found rest counter to 0
        else
            counter +=1 # update counter
        end

        if counter==reset_point
            if best_route_distance < intense_route_distance
                intense_route, intense_route_distance = best_route, best_route_distance
            end
            best_route = shuffle(collect(1:N_cities))
            best_route_distance = calc_distance(best_route)
        end

    end

    if best_route_distance < intense_route_distance
        return best_route , best_route_distance
    else
        return  intense_route, intense_route_distance 
    end
end

##########Testing ############
seed_route = shuffle(collect(1:N_cities))
seed_route_distance = calc_distance(seed_route)
#optimal_distance,optimal_route = three_opt_random(seed_route,seed_route_distance,1000000,coordinates)



optimal_route,optimal_distance = simulated_annealing(coordinates;budget=500000,T=1.0,schedule=0.99,reset_point=1000)

optimal_route = optimal_route .- 1
#######To Do - Try Simulated Annealing #######




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



