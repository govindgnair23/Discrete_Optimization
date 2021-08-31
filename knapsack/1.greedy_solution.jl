# Read in Data
f = open("ks_19_0","r")
input_data = read(f,String)
close(f)

lines = split(input_data,'\n')

items = []

firstLine = split(lines[1])
item_count = parse(Int64,firstLine[1])
capacity = parse(Int64,firstLine[2])

for i ∈ range(1,stop=item_count+1)
    line = lines[i]
    parts = split(line)
    push!(items,(index=i-1,value=parse(Int64,parts[1]),weight = parse(Int64,parts[2])))

end


function greedy(items,capacity)




end