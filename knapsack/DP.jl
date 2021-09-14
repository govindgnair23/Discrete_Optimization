input_file = ARGS[1]

f = open(input_file)
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

using OffsetArrays


#Filter only items that are lighter than the knapsack capacity
items = filter(item-> item.weight <= capacity,items)

"""
Input is list of named tuples and capacity of Knapsack
"""
function DP(items,capacity)
    n_items = length(items)
    #Create a table to hold results of DP
    T = OffsetArray(zeros(Int32,capacity+1,n_items+1),0:capacity,0:n_items)
    #Get minimum weight of items to reduce size of DP table
    #min_weight = minimum(x->x.weight,items)

    for j ∈ range(1,stop=n_items)
        for i ∈ range(1,stop=capacity)

            T[i,j] = T[i,j-1] # Default option is to retain previous item
            if items[j].weight <= i # if  new item is lighter or equal to the current capacity being considered
                option1 = T[i,j-1] # Use previous item 

                option2_1 =  items[j].value    #Select current itenm
                residual_weight = i - items[j].weight    #calculate residual weight after selecting current item           
                option2_2 = residual_weight > 0 ? T[residual_weight,j-1] : 0  # Select best item for residual weight as per DP table         
                T[i,j] = max(option1,option2_1+option2_2)
            end

        end
    end
    # Trace back through T to find selected itens

    i_ = capacity
    j_ = n_items

    selection = []

    while j_ ≠ 0
        if T[i_,j_] == T[i_,j_-1]
            j_ -= 1
        else
            push!(selection,j_)
            i_ = i_ - items[j_].weight
            j_ -= 1
        end
    end

    return selection,T[capacity,n_items] #Last element of the matrix gives value

end
##################Sample Data used for Testing#############################

# items2 = [
#     (index=1,value = 5,weight = 4),
#     (index=2,value = 6,weight = 5),
#     (index=3,value = 3,weight = 2)
# ]

# capacity2 = 9

# selected_items,value = DP(items2,capacity2)

#############################################################################
"""
Function to render output in correct form
item_count: No of items to select from
selected_items: List of selected items by index
value: Value of solution
optimality_flag = 0
"""
function render_output(item_count,selected_items,value,optimality_flag)
    results = fill("0",item_count)
    results[selected_items] .= "1"
    println("$value $optimality_flag")
    println(join(results," "))

end


selected_items,value = DP(items,capacity)

#Get selected items
sel_items = items[selected_items]

#Get indices of selected items
sel_items_i = [item.index for item in sel_items]

render_output(item_count,sel_items_i,value,0)


