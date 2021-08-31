using DataStructures



#Sample problem

items = [
    (index=1,value = 45,weight = 5),
    (index=2,value = 48,weight = 8),
    (index=3,value = 35,weight = 3)
]

K = 10

global N_items = length(items)

global item_weights = map(x-> x.weight, items)
global item_values = map(x-> x.value, items)

#Get cumulative measures for   items
item_cum_weights = cumsum(item_weights)
item_cum_values = cumsum(item_values)

#Sort items by density
items_sorted = sort(items , by = x -> x.value/x.weight,rev=true)

#Get weights of items
item_sorted_weights = map(x-> x.weight, items_sorted)
item_sorted_values = map(x-> x.value, items_sorted)

#Get cumulative measures for  sorted items
item_cum_sorted_weights = cumsum(item_sorted_weights)
item_cum_sorted_values = cumsum(item_sorted_values)

#Get first index where knap sack capacity is exceeded
I = findfirst(x-> x > K, item_cum_sorted_weights)

#Get optimistic estimate -> Convert this to a function
optimistic_value  = item_cum_sorted_values[I-1] #Take all of these items
fraction = (K - item_cum_sorted_weights[I-1])/item_sorted_weights[I]
fractional_value = fraction * item_sorted_values[I]
optimistic_value += fractional_value 

#Get optimistic Knapsack
optimistic_knapsack = push!(fill(1.0,I-1),fraction)

#Function to build an optimistic knapsack given constraints on selecting items

function get_optimistic_estimate(items,K,constraint )
    n_items = length(items)
    dropped = findall(x->x==0,constraint) 
    selected = setdiff(collect(1:n_items),dropped) 

    selected_items = items[selected]
    #Sort items by density
    items_sorted = sort(selected_items , by = x -> x.value/x.weight,rev=true)

    #Get weights of items
    item_sorted_weights = map(x-> x.weight, items_sorted)
    item_sorted_values = map(x-> x.value, items_sorted)

    #Get cumulative measures for  sorted items
    item_cum_sorted_weights = cumsum(item_sorted_weights)
    item_cum_sorted_values = cumsum(item_sorted_values)

   

    if item_cum_sorted_weights[end] < K #if selected items are less than knapsack capacity
        I = length(selected_items)
        optimistic_value = item_cum_sorted_values[I] #select all items
    else
        I = findfirst(x-> x > K, item_cum_sorted_weights)
        optimistic_value  = item_cum_sorted_values[I-1] #Take all of these items
        fraction = (K - item_cum_sorted_weights[I-1])/item_sorted_weights[I]
        fractional_value = fraction * item_sorted_values[I]
        optimistic_value += fractional_value 
    end 

    return optimistic_value
end


##Using Depth First Search
parent = Int[]
#parent = nothing


#Dictionary to hold value,room and estimates
D = Dict(parent =>(value=0,room=K,estimate=optimal_value))

function produce_children(parent::Array{Any,1})
    return [push!(copy(parent),1),push!(copy(parent),0)]
end

#produce_children = (parent::Array{Any,1}) -> [push!(copy(parent),1),push!(copy(parent),0)]

##Function to get updated value,room and estimate for a child produced by produce_children
function evaluate_state(state::Array{Any,1})

    #Get indices corresponding to 1 in the state
    selected = findall(x-> x==1,state)
    dropped = findall(x-> x==0,state)
    

    S_weights = sum(item_weights[selected]) 
    if S_weights < K
        return (value = sum(item_values[selected]),room = K - S_weights,
                    estimate= get_optimistic_estimate(items,K,state ) )
    else
        return (value = nothing, room = K - S_weights,estimate=nothing )
    end

end



######### Development ############

s = Stack{Array{Int64,1}}()
child1,child2= produce_children(parent)
push!(s,child1)
push!(s,child2)

### Convert to Struct for ease of use
struct State    
    state::Array{Any,1} #Also gives the path to the node
    value::Int
    room ::Int
    estimate::Float
end

function produce_children_(parent::State)
    child1,child2 = produce_children(parent.state)
    c1 = evaluate_state(child1)
    c2 = evaluate_state(child2)
    Child1 = State(child1,c1.value,c1.room,c1.estimate)
    Child2 = State(child2,c2.value,c2.room,c2.estimate)
    return Child1,Child2
end



"""
Function to implement depth first search
"""
function DFS2(parent::State)
    explored = Deque{State}()
    to_explore = Deque{State}()
    push!(to_explore,parent)
    
    max_value = 0
    best_solution = nothing
    
    while !isempty(to_explore)
        child = popfirst!(to_explore)

        if (child.value === nothing) 
            continue
        
        elseif  (child.estimate < max_value)
            continue

        elseif length(child.state) == N_items
            if child.value > max_value
                best_solution = child
                max_value = child.value
            end

        else
            child1_,child2_ = produce_children_(child)
            pushfirst!(to_explore,child2_)
            pushfirst!(to_explore,child1_)
        end

        push!(explored,child)
    
    end


    return best_solution


end



##Testing

parent = []

child1,child2 = produce_children(parent)
evaluate_state(child1)
evaluate_state(child2)

child3,child4 = produce_children(child1)
evaluate_state(child3)
evaluate_state(child4)

child5,child6 = produce_children(child4)
evaluate_state(child5)
evaluate_state(child6)

#Test 2

parent = State([],0,0,0.0)
child1,child2 = produce_children_(parent)
child3,child4 = produce_children_(child1)

DFS2(parent)

#Multiple dipstach

function add(x::Int,y::Int)
    return x+y
end

function add(x::String,y::String)
    return x*y
end