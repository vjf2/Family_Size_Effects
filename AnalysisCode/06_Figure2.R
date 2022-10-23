#############
#### Figure 2
#############

library(SocGen)
library(igraph)
library(raster)
library(rgdal)
library(viridisLite)
library(GISTools)

indiv_covars <- read.csv("SharedData/individual_covariates.csv")
focaldata <- read.csv("IntermediateData/affiliations.csv")
coast_polygon <- rgdal::readOGR("SharedData/coastpolygon2020", "coastpolygon2020")

focaldata$residweight <- focaldata$realHWI - focaldata$quantile975
focaldata$residweight <- ifelse(focaldata$residweight < 0, 0, focaldata$residweight)

focaldata$rweight <- focaldata$affiliation975 * focaldata$residweight

affilmat <- dat2mat(focaldata[, c("ID1", "ID2", "rweight")], diag = TRUE)

affilmat[is.na(affilmat)] <- 0

g <- graph.adjacency(affilmat, mode = "undirected", weighted = TRUE, diag = FALSE)

V(g)$degree <- degree(g)

V(g)$vsize <- as.vector(scale(V(g)$degree,
                            center = min(V(g)$degree),
                            scale = diff(range(V(g)$degree))))

cpts <- indiv_covars[, c("centroid_long", "centroid_lat")]
rownames(cpts) <- indiv_covars$Dolphin.ID

cpts <- cpts[which(rownames(cpts) %in% V(g)$name), ]
cpts <- as.matrix(cpts)
cpts <- cpts[V(g)$name, ] # order

# set extent
ylim = c(7138794, 7160609)
xlim = c(760813, 782022)

n = length(V(g))

# add colors and symbols for foraging type and provisioning status

V(g)$fortype <- indiv_covars$ForagingType[match(V(g)$name, indiv_covars$Dolphin.ID)]

V(g)$forcolor[V(g)$fortype == 1] <- viridis(4)[1]
V(g)$forcolor[V(g)$fortype == 2] <- viridis(4)[2]
V(g)$forcolor[V(g)$fortype == 3] <- viridis(4)[4]
V(g)$forcolor[V(g)$fortype == 4] <- viridis(4)[3]
V(g)$forcolor[is.na(V(g)$fortype)] <- "white"

# add triangle shape for provisioned
mytriangle <- function(coords, v = NULL, params) {
        vertex.color <- params("vertex", "color")
        if (length(vertex.color) != 1 && !is.null(v)) {
                vertex.color <- vertex.color[v]
        }
        vertex.size <- 1 / 200 * params("vertex", "size")
        if (length(vertex.size) != 1 && !is.null(v)) {
                vertex.size <- vertex.size[v]
        }

        symbols(x = coords[, 1], y = coords[, 2], bg = vertex.color,
                stars = cbind(vertex.size, vertex.size, vertex.size),
                add = TRUE, inches = FALSE)
}

add_shape("triangle", clip = shapes("circle")$clip,
          plot = mytriangle)

V(g)$status <- indiv_covars$provisioned[match(V(g)$name, indiv_covars$Dolphin.ID)]
V(g)$shape1 <- ifelse(V(g)$status == "Y", "none", "circle")
V(g)$shape2 <- ifelse(V(g)$status == "Y", "triangle", "none")

windows()
# pdf(file="TablesFigures/Figure2.pdf")
par(mar = c(0, 0, 0, 0))
plot(coast_polygon,
     xlim = xlim,
     ylim = ylim,
     col = rgb(0, 0, 0, 50, maxColorValue = 255), border = NA)

# plot edges
plot(g,
     vertex.size = (V(g)$vsize[1:(n)] * 60000 + 17000),
     add = TRUE,
     rescale = FALSE,
     ylim = ylim,
     xlim = xlim,
     vertex.shape = "circle",
     vertex.color = "white",
     edge.arrow.size = 0,
     vertex.label = NA,
     edge.color = adjustcolor("black", alpha.f = 0.4),
     layout = cpts,
     edge.width = edge_attr(g)$weight * 20,
     edge.curved = rep(-.4, length(edge_attr(g)$weight))
)

# make nodes transparent to each other but not to edges
plot(g,
     vertex.size = (V(g)$vsize[1:(n)] * 60000 + 17000),
     add = TRUE,
     rescale = FALSE,
     ylim = ylim,
     xlim = xlim,
     vertex.shape = V(g)$shape1,
     vertex.color = adjustcolor(V(g)$forcolor, alpha.f = 0.7),
     edge.arrow.size = 0,
     vertex.label = NA,
     edge.color = NA,
     layout = cpts,
     edge.width = 0,
     edge.curved = rep(-.4, length(edge_attr(g)$weight))
)

plot(g,
     vertex.size = (V(g)$vsize[1:(n)] * 60000 + 17000),
     add = TRUE,
     rescale = FALSE,
     ylim = ylim,
     xlim = xlim,
     vertex.shape = V(g)$shape2,
     vertex.color = adjustcolor(V(g)$forcolor, alpha.f = 0.7),
     edge.arrow.size = 0,
     vertex.label = NA,
     edge.color = NA,
     layout = cpts,
     edge.width = 0,
     edge.curved = rep(-.4, length(edge_attr(g)$weight))
)

points(772800, 7144424, col = "black", pch = 15, cex = 1.5)

scalebar(10000, xy = c(760348, 7139289), label = c(0, 5, 10), type = "bar", divs = 2, below = "kilometers", font = 2, col = "black", adj = c(0.5, -2))
GISTools::north.arrow(xb = 771348, yb = 7139289, len = 200, lab = "N")

legend("topright", title = "Foraging Category",
       legend = c("Tdpdfor",
                  "Mixed",
                  "Sponging",
                  "Seagrass",
                  "Not assigned",
                  "",
                  "Provisioned"),
       pch = c(rep(21, 5), NA, 24),
       pt.bg = c(viridis(4)[c(2, 1, 4, 3)], "white", NA, "white"),
       pt.cex = c(rep(2, 5), NA, 1.5), bty = "n"
       )

dev.off()
