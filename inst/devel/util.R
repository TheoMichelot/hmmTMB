
# utility for distnames 
update_distnames <- function(name) {
  load("inst/distnames.RData")
  code <- max(distnames$code) + 1
  distnames <- rbind(distnames, data.frame(name = name, code = code))
  save(distnames, file = "inst/distnames.RData")
}

# update documentation
# (this requires temporarily copying distnames.RData to the root of the package)
update_doc <- function() {
    file.copy(from = "inst/distnames.RData", to = "distnames.RData")
    devtools::document()
    file.remove("distnames.RData")
}
