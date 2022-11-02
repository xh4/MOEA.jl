function test_image(img_width, img_height)
    rand(GLubyte, 4, img_width, img_height)
end

function load_image()
    img = load("/mnt/c/Users/XH/Desktop/AwesomeFace.png")'
    img = coloralpha.(img)
    view = channelview(img)
    reinterpret(UInt8, view)
end

function plot_image(fig)
    io = IOBuffer()
    show(io, MIME("image/png"), fig)
    load(io)'
end

function image_uint8(img)
    img = coloralpha.(img)
    view = channelview(img)
    reinterpret(UInt8, view)
end

function fit_image(img, size)
    size = Int(round(size))
    width, height = Base.size(img)
    if width > size
        new_width = size
        new_height = Int(round(size/width*height))
        img = imresize(img, (new_width, new_height))
    end
    img
end

function get_figure_size()
    # CImGui.GetWindowWidth() / 5 * 4
    425
end

g_ImageTexture = Dict{Int,GLuint}()

function CreateImageTexture(image_width, image_height; format=GL_RGBA, type=GL_UNSIGNED_BYTE)
    id = GLuint(0)
    @c glGenTextures(1, &id)
    glBindTexture(GL_TEXTURE_2D, id)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
    glPixelStorei(GL_UNPACK_ROW_LENGTH, 0)
    glTexImage2D(GL_TEXTURE_2D, 0, format, GLsizei(image_width), GLsizei(image_height), 0, format, type, C_NULL)
    g_ImageTexture[id] = id
    return Int(id)
end

function UpdateImageTexture(id, image_data, image_width, image_height; format=GL_RGBA, type=GL_UNSIGNED_BYTE)
    glBindTexture(GL_TEXTURE_2D, g_ImageTexture[id])
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, GLsizei(image_width), GLsizei(image_height), format, type, image_data)
end

function DestroyImageTexture(id)
    id = g_ImageTexture[id]
    @c glDeleteTextures(1, &id)
    delete!(g_ImageTexture, id)
    return true
end

function DestroyImages(ctx)
    for (k,v) in ctx.ImageTexture
        DestroyImageTexture(v)
    end
    return true
end
