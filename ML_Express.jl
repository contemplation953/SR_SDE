import SymbolicRegression: SRRegressor
import MLJ: machine, fit!, predict, report
using MAT

function trans2vect(mat)
    res=Vector{Float64}(undef,length(mat))
    for i in 1:length(mat)
        res[i]=mat[i]
    end
    return res
end

data = matread("D:/ChenXiYuan/workspace/jupyter/Stochastic_ML/H5DOF/data/ML_init_params.mat")

h_val=data["fh_val"][1:15:end]
d_val=data["fd_val"][1:15:end]
mh_val=data["fmh_val"][1:15:end]
sh_val=data["fsh_val"][1:15:end]

X = (y = trans2vect(h_val), x = trans2vect(d_val))
y = trans2vect(mh_val)
z = trans2vect(sh_val)

model = SRRegressor(
    niterations=50,
    binary_operators=[+, -, *],
    # batching=true,
    # batch_size=10000
)

mach1 = machine(model, X, y)
@time fit!(mach1)

mach2 = machine(model, X, z)
@time fit!(mach2)
