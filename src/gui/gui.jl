glsl_version = 130

function __init__()
    @static if Sys.isapple()
        # OpenGL 3.2 + GLSL 150
        global glsl_version = 150
        GLFW.WindowHint(GLFW.CONTEXT_VERSION_MAJOR, 3)
        GLFW.WindowHint(GLFW.CONTEXT_VERSION_MINOR, 2)
        GLFW.WindowHint(GLFW.OPENGL_PROFILE, GLFW.OPENGL_CORE_PROFILE) # 3.2+ only
        # GLFW.WindowHint(GLFW.OPENGL_FORWARD_COMPAT, GLFW.GL_TRUE) # required on Mac
    else
        # OpenGL 3.0 + GLSL 130
        global glsl_version = 130
        GLFW.WindowHint(GLFW.CONTEXT_VERSION_MAJOR, 3)
        GLFW.WindowHint(GLFW.CONTEXT_VERSION_MINOR, 0)
        # GLFW.WindowHint(GLFW.OPENGL_PROFILE, GLFW.OPENGL_CORE_PROFILE) # 3.2+ only
        # GLFW.WindowHint(GLFW.OPENGL_FORWARD_COMPAT, GL_TRUE) # 3.0+ only
    end
end

error_callback(err::GLFW.GLFWError) = @error "GLFW ERROR: code $(err.code) msg: $(err.description)"

function init_renderer(width, height, title::AbstractString)
    # setup GLFW error callback
    GLFW.SetErrorCallback(error_callback)

    # create window
    window = GLFW.CreateWindow(width, height, title)
    @assert window != C_NULL
    GLFW.MakeContextCurrent(window)
    GLFW.SwapInterval(1)  # enable vsync

    # setup Dear ImGui context
    ctx = CImGui.CreateContext()

    # setup Dear ImGui style
    CImGui.StyleColorsDark()
    # CImGui.StyleColorsClassic()
    # CImGui.StyleColorsLight()

    # setup Platform/Renderer bindings
    CImGui.ImGui_ImplGlfw_InitForOpenGL(window, true)
    CImGui.ImGui_ImplOpenGL3_Init(glsl_version)

    return window, ctx
end

gui_task = nothing

function renderloop(window, ctx, ui=()->nothing, hotloading=false)
    try
        while !GLFW.WindowShouldClose(window)
            GLFW.PollEvents()
            CImGui.ImGui_ImplOpenGL3_NewFrame()
            CImGui.ImGui_ImplGlfw_NewFrame()
            CImGui.NewFrame()

            hotloading ? Base.invokelatest(ui) : ui()

            CImGui.Render()
            GLFW.MakeContextCurrent(window)
            display_w, display_h = GLFW.GetFramebufferSize(window)
            glViewport(0, 0, display_w, display_h)
            glClearColor(0.2, 0.2, 0.2, 1)
            glClear(GL_COLOR_BUFFER_BIT)
            CImGui.ImGui_ImplOpenGL3_RenderDrawData(CImGui.GetDrawData())

            GLFW.MakeContextCurrent(window)
            GLFW.SwapBuffers(window)
            yield()
        end
    catch e
        @error "Error in renderloop!" exception=e
        Base.show_backtrace(stderr, catch_backtrace())
    finally
        CImGui.ImGui_ImplOpenGL3_Shutdown()
        CImGui.ImGui_ImplGlfw_Shutdown()
        CImGui.DestroyContext(ctx)
        GLFW.DestroyWindow(window)
    end
end

function render(ui; width=1280, height=720, title::AbstractString="Demo", hotloading=false)
    window, ctx = init_renderer(width, height, title)
    GC.@preserve window ctx begin
        t = @async renderloop(window, ctx, ui, hotloading)
    end
    return t
end

function gui()
    show_playground_window = true
    show_experiment_window = true
    render_playground_window = playground_window()
    render_experiment_window = experiment_window()
    if gui_task !== nothing
        # ex = InterruptException()
        # Base.throwto(gui_task, ex)
    end
    Plots.theme(:dark)
    global gui_task = Threads.@spawn Threads.@sync render(width = 600, height = 800, title = "MOEA") do
        CImGui.Begin("Start")
        if CImGui.Button("Playground")
            show_playground_window = !show_playground_window
        end
        if show_playground_window
            render_playground_window()
        end
        CImGui.SameLine()
        if CImGui.Button("Experiment")
            show_experiment_window = !show_experiment_window
        end
        if show_experiment_window
            render_experiment_window()
        end
        CImGui.End()
    end
end
