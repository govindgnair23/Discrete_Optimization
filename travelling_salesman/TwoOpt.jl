# Implement the 2-Opt technique with simulated annealing
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

#Function to calculate distance given coordinates
distance(city1,city2) = sqrt((city1[1]-city2[1])^2 + (city1[2]-city2[2])^2)

#Function to calculate length of two edges in a sequence
#Computes distance node1-node2 and node3-node4
function distance2opt((node1,node2,node3,node4),)
     return distance(coordinates[node1],coordinates[node2]) + distance(coordinates[node3],coordinates[node4])
end

"""
Function to do local search with 2-Opt
"""
function do_two_opt(sequence::Array{Int64,1})
    L = length(sequence)
    #Select two edges by selecting 2 elements from the first n-1 elements of the sequence
    edge1_node1 = sample(1:length(sequence),1)[1]
    diff = edge1_node1 == 1 ?  [edge1_node1,edge1_node1+1,L] : [edge1_node1-1,edge1_node1,edge1_node1+1]
    edge2_node1 = sample(setdiff(collect(1:length(sequence)),diff),1)[1]
    edge1_node2 = sequence[(edge1_node1+1)%L]
    edge2_node2 = sequence[(edge2_node1+1)%L]

    #Get original distance between selected edges
    original_distance = distance2opt((edge1_node1,edge1_node2,edge2_node1,edge2_node2))

    possible_edges = [(edge1_node1,edge2_node1,edge1_node2,edge2_node2),(edge1_node1,edge2_node2,edge1_node2,edge2_node1)]
    #Evaluate edges and select best edge
    #Evaluate option1
    option1 = distance2opt(possible_edges[1])
    ##valuate option 2
    option1 = distance2opt(possible_edges[2])








    return sequence
end

