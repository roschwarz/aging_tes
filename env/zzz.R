.onLoad <- function(libname, pkgname) {
    message("Loading project environment: example")
    
    # New isolated Env
    env <- new.env(parent = emptyenv())
    
    # Load modules directly
    sys.source(system.file("R/setup.R", package = pkgname), envir = env)
    sys.source(system.file("R/constants.R", package = pkgname), envir = env)
    
    # export to global environment
    list2env(as.list(env), envir = parent.env(environment()))
}