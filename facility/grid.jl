#Assign each customer to nearest facility

using LinearAlgebra
using Plots
import HiGHS
using JuMP

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
N_facilties = parse(Int64,firstLine[1])
N_customers = parse(Int64,firstLine[2])

for i ∈ range(2,stop=N_facilties+1)
    line = lines[i]
    parts = split(line)
    push!(f,parse(Float64,parts[1]))
    push!(q,parse(Int64,parts[2]))
    push!(Xf,parse(Float64,parts[3]))
    push!(Yf,parse(Float64,parts[4]))
end

for i ∈ range(N_facilties+2,stop=N_facilties+N_customers+1)
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
global C = zeros(N_customers,N_facilties) # Computed for all customers and facilities
for i in 1:N_customers
    for j in 1:N_facilties
        C[i, j] = LinearAlgebra.norm([Xc[i] - Xf[j], Yc[i] - Yf[j]], 2)
    end
end

#Find closest facility to each customer

assignment  = map(eachrow(C)) do r
    argmin(r)
end

y = zero(assignment)
y[assignment] .= 1


#Get max coordinates to consider
X_max = ceil(maximum([maximum(Xc),maximum(Xf)]))
Y_max = ceil(maximum([maximum(Yc),maximum(Yf)]))


N = 5 #No of cuts- creates N-1 partitions

X_cuts = collect(range(0,X_max,length=5))
Y_cuts = collect(range(0,Y_max,length=5))


###Loop through all grids and solve

for i ∈ 1:N-1
    for j ∈ 1:N-1
        X_grid = findall( x -> (X_cuts[i]<=x<=X_cuts[i+1]),Xc)
        Y_grid = findall( y -> (Y_cuts[j]<=y<=Y_cuts[j+1]),Yc)
        grid_customers = intersect(X_grid,Y_grid)

        assignment = assign_facilities(grid_customers,assignment,C)
    end
end

##Function to optimally assign a set of customers
"""
customers:: Set of customers to be assigned (1 D Array)
assignment:: current assignment mapping each customer to a facility(1 D array)
c:: cost matrix
"""
function assign_facilities(customers::Vector{Int64},assignment::Vector{Int64})

    #Get facilities assigned to these customers as per current solution
    facilities = assignment[customers]

    #Get all customers assigned to these facilties
    custs = findall(i -> (i ∈ facilities),assignment)

    #Add these customers to the customers in the grid
    custs = unique(union(customers,custs))

    n_c = length(custs) # Number of customers in this sub-problem
    n_f = length(facilities) #Number of facilities in this sub-problem
    #Cost matrix for just these customers and facilities
    c_ = C[custs,facilities]

    a_ = a[custs]
    q_ = q[facilities]
    f_ = f[facilities]

    cfl = Model(HiGHS.Optimizer)

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
function calc_cost(assignment::Vector{Int64},N_facilties)
    #Get all opened facilities
    opened_facilities = unique(assignment)
    #Convert assignment to a matrix to get transportation costs with simple multiplication
    x = M.OneHotMat(assignment,N_facilties)
    #Get transoration costs
    TC = sum(C.*x) 
    #Fixed Costs
    FC = sum(f[opened_facilities])
    return TC+FC
end

#######Try solving MIP with just customers and facilties in just one grid
#Get all customers in the first grid

X_grid1 = findall( i -> (X_cuts[1]<=i<=X_cuts[2]),Xc)
Y_grid1 = findall( i -> (Y_cuts[1]<=i<=Y_cuts[2]),Yc)
grid1_customers = intersect(X_grid1,Y_grid1)

#Get facilities assigned to these customers
facilities = assignment[grid1_customers]

#Get all customers assigned to these facilties
custs = findall(i -> (i ∈ facilities),assignment)

#Add these customers to the customers in the grid
custs = unique(union(grid1_customers,custs))


#Constuct a cost matrix for just these customers and facilties
n_c = length(custs)
n_f = length(facilities)
c = zeros(n_c,n_f)

for i in 1:n_c
    for j in 1:n_f
        c[i, j] = LinearAlgebra.norm([Xc[custs[i]] - Xf[facilities[j]], Yc[custs[i]] - Yf[facilities[j]]], 2)
    end
end

a_ = a[custs]
q_ = q[facilities]
f_ = f[facilities]

cfl = Model(HiGHS.Optimizer)

@variable(cfl, y[1:n_f],Bin);
@variable(cfl, x[1:n_c,1:n_f],Bin);

#Each client is served from only one warehouse
@constraint(cfl, client_service[i in 1:n_c], sum(x[i,:]) == 1);

#Capacity constraint
@constraint(cfl, capacity, x'a_ .<= (q_.*y));

#Objective
@objective(cfl,Min, f_'y + sum(c .* x));


optimize!(cfl)

obj_value = objective_value(cfl)


#Visualize these customers and facilties
Xc1 = Xc[custs]
Yc1 = Yc[custs]
Xf1 = Xf[facilities]
Yf1 = Yf[facilities]

#Update original assignment for these customers
x_ = value.(x) .> 1 - 1e-5
#Get index of column assigned 1
i_ = map(i -> i[2],argmax(x_,dims=2))

#Test function
assign_facilities(grid1_customers,assignment,C)

#Updates assignments to these customers
assignment[custs] .= facilities[i_[:]]

plot_solution(x,y,Xc1,Yc1,Xf1,Yf1)








plot_locations(Xc1,Yc1,Xf1,Yf1)



function plot_locations(Xc,Yc,Xf,Yf)
    p = Plots.scatter(
        Xc,
        Yc,
        label = nothing,
        markershape = :circle,
        markercolor = :blue,
    )

    Plots.scatter!(
    Xf,
    Yf,
    label = nothing,
    markershape = :rect,
    markerstrokecolor = :red,
    markerstrokewidth = 2,
    )

    p
end

function plot_solution(x,y,Xc,Yc,Xf,Yf)
    x_ = value.(x) .> 1 - 1e-5;
    y_ = value.(y) .> 1 - 1e-5;
    N_customers,N_facilities = size(x) 

    p = Plots.scatter(
        Xc,
        Yc,
        label = nothing,
        markershape = :circle,
        markercolor = :blue,
    )


    #mc = [(y_[j] ? :red : :white) for j in 1:N_facilities]



    Plots.scatter!(
        Xf[y_],
        Yf[y_],
        label = nothing,
        markershape = :rect,
        markercolor = :red,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )

    for i in 1:N_customers
        for j in 1:N_facilities
            if x_[i,j] == 1
                Plots.plot!(
                    [Xc[i],Xf[j]],
                    [Yc[i], Yf[j]],
                    color = :black,
                    label = nothing,
                )
                break
            end
        end
    end

    p

end