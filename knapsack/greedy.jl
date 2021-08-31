
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



    
"""
Input is list of named tuples and capacity of Knapsack
"""
function greedy(items,capacity)
    value_density = [item.value/item.weight for item in items]
    value_rank = sortperm(value_density,rev=true)

    C =  capacity
    K = [] #Knapsack
    V = 0 #value

    for i ∈ value_rank
        i_weight = items[i].weight
        if i_weight <= C
            push!(K,i)
            C -= i_weight
            V += items[i].value 
        end
    end
    return K,V
end

selected_items,value = greedy(items,capacity)


function render_output(item_count,selected_items,value,optimality_flag)
    results = fill("0",item_count)
    results[selected_items] .= "1"
    println("$value $optimality_flag")
    println(join(results," "))

end

render_output(item_count,selected_items,value,0)

