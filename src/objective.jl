mutable struct Objective{TF,TX,TY} <: AbstractObjective
    # 函数
    f::TF
    # 输入
    x::TX
    # 输出
    y::TY
    # 调用次数
    f_calls::Int
end

"""
    Objective(f, x[, y])
    
f: 目标函数
x: 输入值
y: 输出值
"""
function Objective(
    f::TF,
    x::AbstractArray,
    y::Union{Real,AbstractArray{<:Real}} = zero(f(x))
) where {TF}
    Objective{TF,typeof(x),typeof(y)}(f, x, y, 0)
end

function Objective(
    f::TF,
    x::Individual,
    y::Union{Real,AbstractArray{<:Real}} = zero(f(variable(x)))
) where {TF}
    Objective{TF,typeof(variable(x)),typeof(y)}(f, variable(x), y, 0)
end

f_calls(obj::Objective) = obj.f_calls

value(obj::Objective) = obj.y

"""
    ismultiobjective(objfun)

判断是否为多目标
"""
ismultiobjective(obj) = obj.y isa AbstractArray

function value(obj::Objective, x)
    obj.f_calls += 1
    obj.f(x)
end

function value(obj::Objective, x::Individual)
    value(obj, variable(x))
end

function value!(obj::Objective, x)
    obj.y = value(obj, x)
    obj.y
end

function value!(obj::Objective, x::Individual)
    obj.y = value(obj, x)
    x.objective = obj.y
    obj.y
end

function value!(
    obj::Objective,
    xs::Union{AbstractVector, Population},
)
    map((x) -> value!(obj, x), xs)
end
