#Assign each customer to nearest facility


file = open(input_file)
input_data = read(file,String)
close(file)

lines = split(input_data,'\n')

#Fixed costs
f = []
#Capacities
q = []

#Demands from each of M customers
a = []

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
c = zeros(N_customers,N_facilties)
for i in 1:N_customers
    for j in 1:N_facilties
        c[i, j] = LinearAlgebra.norm([Xc[i] - Xf[j], Yc[i] - Yf[j]], 2)
    end
end

#Find closest facility to each customer

assignment  = map(eachrow(c)) do r
    argmin(r)
end

y = zero(assignment)
y[assignment] .= 1

"""
Function to plot assignment of warehouses to customers
Args:
    assignment: Vector assigning customers to a facility.
"""
function plot_solution(assignment,y,Xc,Yc,Xf,Yf)
    y_ = value.(y) .> 1 - 1e-5;

    p = Plots.scatter(
        Xc,
        Yc,
        label = nothing,
        markershape = :circle,
        markercolor = :blue,
    )


    mc = [(y_[j] ? :red : :white) for j in 1:N_facilties]



    Plots.scatter!(
        Xf,
        Yf,
        label = nothing,
        markershape = :rect,
        markercolor = mc,
        markerstrokecolor = :red,
        markerstrokewidth = 2,
    )

    for i in 1:N_customers
        assigned_facility = assignment[i]
                Plots.plot!(
                    [Xc[i],Xf[assigned_facility]],
                    [Yc[i], Yf[assigned_facility]],
                    color = :black,
                    label = nothing,
                )
    end

    p

end

plot_solution(assignment,y,Xc,Yc,Xf,Yf)