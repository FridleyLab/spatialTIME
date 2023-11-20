library(hexSticker)
library(magick)

# sticker ---------------------------------------------------------------------

# sysfonts::font_add_google("Fjalla One")

img <- image_read("man/hex_sticker/tma_figure.png")

sticker(img, package="spatialTIME",
        p_size=6, p_color = "#40196A", p_y = 1,
        s_x=1, s_y=1, s_width=1.7, s_height=1.7,
        h_fill="white", h_color="#40196A", 
        url = "https://github.com/FridleyLab/spatialTIME", u_color = "#40196A", u_size = 1.1,
        filename="man/hex_sticker/hex.png")