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
    CImGui.GetWindowWidth() / 5 * 4
end
