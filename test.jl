using Test
using Gauss

include("src/bsp.jl")

gauss!(lgs3a)
rückwärtsauflösen!(lgs3a)

gauss!(lgs3b)
gerade(lgs3b)

gauss!(lgs3c)
gerade(lgs3b, var = "t")

gauss!(lgs10a)
gerade(lgs10a)

gauss!(lgs10b)
rückwärtsauflösen!(lgs10b)

gauss!(lgs10c)
gerade(lgs10c)

gauss!(lgs11a)
rückwärtsauflösen!(lgs11a)

gauss!(lgs11b)

gauss!(lgs11c)
