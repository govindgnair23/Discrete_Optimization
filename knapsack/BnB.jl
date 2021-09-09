input_data = ARGS[1]

f = open(input_data)
input_data = read(f,String)
close(f)

lines = split(input_data,'\n')

items = []

firstLine = split(lines[1])
item_count = parse(Int64,firstLine[1])
capacity = parse(Int64,firstLine[2])

for i ∈ range(2,stop=item_count+1)
    line = lines[i]
    parts = split(line)
    push!(items,(index=i-1,value=parse(Int64,parts[1]),weight = parse(Int64,parts[2])))

end

# Remove any items with weight greater than capacity of the knapsack
items= filter(item-> item.weight <= capacity,items)



using DataStructures



global N_items = length(items)

global item_weights = map(x-> x.weight, items)
global item_values = map(x-> x.value, items)



#Function to build an optimistic knapsack given constraints on selecting items
#Constraint is ewuivalent to state or how many items have been selected/dropoed
function get_optimistic_estimate(items,capacity,constraint )
    n_items = length(items)
    dropped = findall(x->x==0,constraint) 
    selected = setdiff(collect(1:n_items),dropped) 

    selected_items = items[selected]
    #Sort items by value density
    items_sorted = sort(selected_items , by = x -> x.value/x.weight,rev=true)

    #Get weights of items
    item_sorted_weights = map(x-> x.weight, items_sorted)
    item_sorted_values = map(x-> x.value, items_sorted)

    #Get cumulative measures for  sorted items
    item_cum_sorted_weights = cumsum(item_sorted_weights)
    item_cum_sorted_values = cumsum(item_sorted_values)

   

    if item_cum_sorted_weights[end] < capacity #if selected items are less than knapsack capacity
        I = length(selected_items)
        optimistic_value = item_cum_sorted_values[I] #select all items
    else
        I = findfirst(x-> x > capacity, item_cum_sorted_weights)
        optimistic_value  = item_cum_sorted_values[I-1] #Take all of these items
        fraction = (capacity - item_cum_sorted_weights[I-1])/item_sorted_weights[I]
        fractional_value = fraction * item_sorted_values[I]
        optimistic_value += fractional_value 
    end 

    return optimistic_value
end


#Given a state/parent in array for produce children
function produce_children(parent::Array{Int64,1})
    return [push!(copy(parent),1),push!(copy(parent),0)]
end



##Function to get updated value,room and estimate for a child produced by produce_children
function evaluate_state(state::Array{Int64,1})

    #Get indices corresponding to 1 in the state
    selected = findall(x-> x==1,state)
    dropped = findall(x-> x==0,state)
    

    S_weights = sum(item_weights[selected]) 
    if S_weights < capacity
        return (value = sum(item_values[selected]),room = capacity - S_weights,
                    estimate= get_optimistic_estimate(items,capacity,state ) )
    else
        return (value = nothing, room = capacity - S_weights,estimate=nothing )
    end

end



### Convert to a Struct for ease of use
struct State    
    state::Array{Int64,1} #Also gives the path to the node
    value::Union{Int64,Nothing}
    room ::Int
    estimate::Union{Float64,Nothing}
end


function produce_children(parent::State)
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
            child1_,child2_ = produce_children(child)
            pushfirst!(to_explore,child2_)
            pushfirst!(to_explore,child1_)
        end

        push!(explored,child)
    
    end


    return best_solution


end





"""
Function to render output in correct form
item_count: No of items to select from
selected_items: Boolean vector represneting items selected
value: Value of solution
optimality_flag = 0
"""
function render_output(item_count,selected_items,value,optimality_flag)
    println("$value $optimality_flag")
    println(join(selected_items," "))
end

parent = State(Int64[],0,0,0.0)
solution = DFS2(parent)
render_output(item_count,solution.state,solution.value,0)
