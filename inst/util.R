# utility for distnames 

update_distnames <- function(name) {
  load("inst/distnames.RData")
  code <- max(distnames$code) + 1
  distnames <- rbind(distnames, data.frame(name = name, code = code))
  save(distnames, file = "inst/distnames.RData")
}
