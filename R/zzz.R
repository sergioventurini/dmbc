.onAttach <- function(lib, pkg) {
	packageStartupMessage(sprintf("Package %s (%s) loaded.
To cite, type citation(\"%s\")", pkg, packageDescription(pkg)$Version, pkg))
}

.onUnload <- function(lib) {
	library.dynam.unload("dmbc", lib)
}
