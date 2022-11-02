mutable struct Objective <: AbstractObjective
    # 函数
    f
    # 输入
    x
    # 输出
    y
    # 调用次数
    f_calls::Int
end

"""
    Objective(f, x[, y])
    
f: 目标函数
x: 输入值
y: 输出值
"""
function Objective(f, x, y)
    Objective(f, x, y, 0)
end

function Objective(f, x)
    Objective(f, x, f(x), 0)
end

function Objective(f, x::AbstractIndividual, y)
    Objective(f, variables(x), y, 0)
end

function Objective(f, x::AbstractIndividual)
    Objective(f, variables(x), f(variables(x)), 0)
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

function value(obj::Objective, x::AbstractIndividual)
    value(obj, variables(x))
end

function value!(obj::Objective, x)
    obj.y = value(obj, x)
    obj.y
end

function value!(obj::Objective, x::AbstractIndividual)
    obj.y = value(obj, x)
    x.objectives = obj.y
    obj.y
end

function value!(
    obj::Objective,
    xs::Union{AbstractVector, Population},
)
    map((x) -> value!(obj, x), xs)
end
