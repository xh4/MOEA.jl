
db = SQLite.DB(joinpath(homedir(), "MOEA.sqlite3"))

function init_db()
    global db = SQLite.DB(joinpath(homedir(), "MOEA.sqlite3"))

    if !("settings" in table_names())
        schema = Tables.Schema(["key", "value"], [String, Any])
        SQLite.createtable!(db, "settings", schema, temp=false, ifnotexists=true)
    end

    if !("runs" in table_names())
        schema = Tables.Schema(["problem", "algorithm", "igd", "hv", "time"], [String, String, Number, Number, Number])
        SQLite.createtable!(db, "runs", schema, temp=false, ifnotexists=true)
    end
end

function ensure_db()
    if db.handle == C_NULL
        init_db()
    end
end

function table_names()
    ensure_db()
    df = DBInterface.execute(db, "SELECT * FROM sqlite_master WHERE type = 'table'") |> DataFrame
    df.name
end

function get_setting(key::AbstractString, default_value=nothing)
    ensure_db()
    df = DBInterface.execute(db, "SELECT value FROM settings WHERE key = ?", [key]) |> DataFrame
    if isempty(df.value)
        try
            deserialize(IOBuffer(df.value[1]))
        catch
            DBInterface.execute(db, "DELETE FROM settings WHERE key = ?", [key])
            default_value
        end
    else
        default_value
    end
end

function has_setting(key::AbstractString)
    ensure_db()
    df = DBInterface.execute(db, "SELECT key FROM settings WHERE key = ?", [key]) |> DataFrame
    !isempty(df.key)
end

function set_setting(key::AbstractString, value)
    ensure_db()
    io = IOBuffer()
    serialize(io, value)
    v = take!(io)
    if has_setting(key)
        DBInterface.execute(db, "UPDATE settings SET value = ? WHERE key = ?", [v, key])
    else
        DBInterface.execute(db, "INSERT INTO settings (key, value) VALUES (?, ?)", [key, v])
    end
end

function save_run(problem, algorithm, igd, hv)
    ensure_db()
    DBInterface.execute(db, "INSERT INTO runs (problem, algorithm, igd, hv, time) VALUES (?, ?, ?, ?, ?)", 
                            [identifier(problem), identifier(algorithm), igd, hv, time()])
end

function fetch_runs(problem, algorithm, n = 30)
    ensure_db()
    DBInterface.execute(db, "SELECT igd, hv FROM runs WHERE problem = ? AND algorithm = ? ORDER BY time DESC LIMIT ?", [identifier(problem), identifier(algorithm), n]) |> DataFrame
end

function clear_runs(; algorithm = nothing, problem = nothing)
    stat = "DELETE FROM RUNS"
    if algorithm !== nothing || problem !== nothing
        stat *= " WHERE"
        if algorithm !== nothing
            algorithm = identifier(algorithm)
            stat *= " algorithm = ?"
        end
        if problem !== nothing
            problem = identifier(problem)
            if algorithm !== nothing
                stat *= " AND"
            end
            stat *= " problem = ?"
        end
    end
    vals = filter(x -> !isnothing(x), [algorithm, problem])
    DBInterface.execute(db, stat, vals)
end