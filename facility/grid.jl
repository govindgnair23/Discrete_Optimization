##Solve problem by breaking up entire region into grids and solving sub-problem in each grid.

#Assign each customer to nearest facility
input_file = ARGS[1]


using LinearAlgebra
using Plots
import HiGHS
using JuMP
using Printf

file = open(input_file)
input_data = read(file,String)
close(file)

lines = split(input_data,'\n')

#Fixed costs
global f = []
#Capacities
global q = []

#Demands from each of M customers
global a = []

#X and Y coordinates of customers and facilties
Xc,Yc,Xf,Yf = [],[],[],[]

firstLine = split(lines[1])
#Number of facilites
N_facilities = parse(Int64,firstLine[1])
N_customers = parse(Int64,firstLine[2])

for i ∈ range(2,stop=N_facilities+1)
    line = lines[i]
    parts = split(line)
    push!(f,parse(Float64,parts[1]))
    push!(q,parse(Int64,parts[2]))
    push!(Xf,parse(Float64,parts[3]))
    push!(Yf,parse(Float64,parts[4]))
end

for i ∈ range(N_facilities+2,stop=N_facilities+N_customers+1)
    line = lines[i]
    parts = split(line)
    push!(a,parse(Int64,parts[1]))
    push!(Xc,parse(Float64,parts[2]))
    push!(Yc,parse(Float64,parts[3]))
end

# Convert to appropriate format
f = convert(Array{Float64},f)
q = convert(Array{Int64},q)
a = convert(Array{Int64},a)

#Transportation Costs
global C = zeros(N_customers,N_facilities) # Computed for all customers and facilities
for i in 1:N_customers
    for j in 1:N_facilities
        C[i, j] = LinearAlgebra.norm([Xc[i] - Xf[j], Yc[i] - Yf[j]], 2)
    end
end

#######Greedy feasible solution ###########
#Find closest facility to each customer, if the facility is already full, assign to next closest facility
facility_loads = zeros(N_facilities)
assignment = zeros(Int64,N_customers)

for i ∈ 1:N_customers
    #Get nearest facilities to customer i
    closest_facilities = sortperm(C[i,:])
    assigned = false
    j = 1 # Start by assigning closest facility
    while !assigned #while the customer has not been assigned to a facility
        assigned_facility = closest_facilities[j]  #Start by assigning closest facility
        if a[i] <= q[assigned_facility] - facility_loads[assigned_facility] #If there is enough room to accomodate new customer
            assignment[i] = assigned_facility
            facility_loads[assigned_facility] +=a[i] # Add to facilitiy's assigned facility_loads
            assigned = true
        else
            j += 1 #Try next clostest facility
        end
    end

end



y = zero(1:N_facilities)



#Get max coordinates to consider
X_max = ceil(maximum([maximum(Xc),maximum(Xf)]))
Y_max = ceil(maximum([maximum(Yc),maximum(Yf)]))

if N_customers + N_facilities >= 1000
    N = 20 #No of cuts- creates N-1 partitions
elseif N_customers + N_facilities >= 500
    N =  10
elseif N_customers + N_facilities > 250
    N = 5
else
    N = 2 #whole region will be considered 1 grid
end


X_cuts = collect(range(0,X_max,length=N))
Y_cuts = collect(range(0,Y_max,length=N))

##Function to optimally assign a set of customers
"""
customers:: Set of customers to be assigned (1 D Array)
assignment:: current assignment mapping each customer to a facility(1 D array)
c:: cost matrix
"""
function assign_facilities(customers::Vector{Int64},assignment::Vector{Int64},a,f,q,C)

    #Get facilities assigned to these customers as per current solution
    facilities = unique(assignment[customers])

    #Get all customers assigned to these facilties
    custs = findall(i -> (i ∈ facilities),assignment)

    #Add these customers to the customers in the grid
    custs = unique(union(customers,custs))

    n_c = length(custs) # Number of customers in this sub-problem
    n_f = length(facilities) #Number of facilities in this sub-problem

    if n_f <= 1 # If there is only one facility for a customer to be assigned to there is no optimization ptoblem
        return assignment
    end
    #Cost matrix for just these customers and facilities
    c_ = C[custs,facilities]

    a_ = a[custs]
    q_ = q[facilities]
    f_ = f[facilities]

    cfl = Model(HiGHS.Optimizer)
    set_optimizer_attribute(cfl,"log_to_console",false)

    @variable(cfl, y[1:n_f],Bin);
    @variable(cfl, x[1:n_c,1:n_f],Bin);

    #Each client is served from only one warehouse
    @constraint(cfl, client_service[i in 1:n_c], sum(x[i,:]) == 1);

    #Capacity constraint
    @constraint(cfl, capacity, x'a_ .<= (q_.*y));

    #Objective
    @objective(cfl,Min, f_'y + sum(c_ .* x));


    optimize!(cfl)

    #Update original assignment for these customers
    x_ = value.(x) .> 1 - 1e-5
    #Get index of column assigned 1
    i_ = map(i -> i[2],argmax(x_,dims=2))

    #Updates assignments to these customers
    assignment[custs] .= facilities[i_[:]]

    return assignment
end

module M
    struct OneHotMat{T} <:AbstractMatrix{T} 
        v::AbstractVector{T}
        f::Int64
    end

    function Base.size(M::OneHotMat)
        return length(M.v),M.f
    end

    function Base.getindex(M::OneHotMat,i,j)
        return Int(M.v[i]==j)
    end
end

#Function to calculate total cost from an assignment
function calc_cost(assignment::Vector{Int64},N_facilities,C)
    #Get all opened facilities
    opened_facilities = unique(assignment)
    #Convert assignment to a matrix to get transportation costs with simple multiplication
    x = M.OneHotMat(assignment,N_facilities)
    #Get transoration costs
    TC = sum(C.*x) 
    #Fixed Costs
    FC = sum(f[opened_facilities])
    return TC+FC,x
end



function solve_problem(X_cuts,Y_cuts,assignment,a,f,q)
    N = length(X_cuts)
    ###Loop through all grids and solve
    for i ∈ 1:N-1
        for j ∈ 1:N-1
            X_grid = findall( x -> (X_cuts[i]<=x<=X_cuts[i+1]),Xc)
            Y_grid = findall( y -> (Y_cuts[j]<=y<=Y_cuts[j+1]),Yc)
            grid_customers = intersect(X_grid,Y_grid)
            new_assignment = assign_facilities(grid_customers,assignment,a,f,q,C)
            assignment = new_assignment
        end
    end
    return assignment
end



final_assignment = solve_problem(X_cuts,Y_cuts,assignment,a,f,q)



#Get final assignment in matrix form and objective value
cost,x = calc_cost(final_assignment,N_facilities,C)

#Set only opened facilities to 1
y[unique(final_assignment)] .= 1

#########Test whether any facility constraints were violated ######

# facility_loads2 =  zeros(N_facilities)
# for i ∈ 1:N_facilities
#     custs_i = findall(c -> c== i, final_assignment)
#     facility_loads2[i] = sum(a[custs_i])
# end

# sum(facility_loads2 .> q)


##########################################################################
# function plot_locations(Xc,Yc,Xf,Yf)
#     p = Plots.scatter(
#         Xc,
#         Yc,
#         label = nothing,
#         markershape = :circle,
#         markercolor = :blue,
#     )

#     Plots.scatter!(
#     Xf,
#     Yf,
#     label = nothing,
#     markershape = :rect,
#     markerstrokecolor = :red,
#     markerstrokewidth = 2,
#     )

#     p
# end

# function plot_solution(x,y,Xc,Yc,Xf,Yf)
#     x_ = value.(x) .> 1 - 1e-5;
#     y_ = value.(y) .> 1 - 1e-5;
#     N_customers,N_facilities = size(x) 

#     p = Plots.scatter(
#         Xc,
#         Yc,
#         label = nothing,
#         markershape = :circle,
#         markercolor = :blue,
#     )


#     #mc = [(y_[j] ? :red : :white) for j in 1:N_facilities]



#     Plots.scatter!(
#         Xf[y_],
#         Yf[y_],
#         label = nothing,
#         markershape = :rect,
#         markercolor = :red,
#         markerstrokecolor = :red,
#         markerstrokewidth = 2,
#     )

#     for i in 1:N_customers
#         for j in 1:N_facilities
#             if x_[i,j] == 1
#                 Plots.plot!(
#                     [Xc[i],Xf[j]],
#                     [Yc[i], Yf[j]],
#                     color = :black,
#                     label = nothing,
#                 )
#                 break
#             end
#         end
#     end

#     p

# end


#plot_locations(Xc,Yc,Xf,Yf)
#plot_solution(x,y,Xc,Yc,Xf,Yf)

#############################################################################
"""
Function to render output in correct form
item_count: No of items to select from
selected_items: List of selected items by index
value: Value of solution
optimality_flag = 0
"""
function render_output(assignment,cost;optimality_flag=0)
    @printf("%0.2f ",cost)
    println("$optimality_flag")
    println(join(assignment," "))
end

final_assignment = final_assignment .- 1 # Solution requires indexing from 0
render_output(final_assignment,cost)

