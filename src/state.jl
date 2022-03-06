abstract type AbstractOptimizerState end

"""
    value(state)

返回当前状态的极小值
"""
value(state::AbstractOptimizerState) = error("`value` is not implemented for $(state).")

"""
    minimizer(state)

返回当前状态的决策变量
"""
minimizer(state::AbstractOptimizerState) = error("`minimizer` is not implemented for $(state).")

"""
    terminate(state)

判断是否提前终止
"""
terminate(state::AbstractOptimizerState) = false