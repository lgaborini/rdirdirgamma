
# Detach the DLL when unloading the library
.onUnload <- function(libpath) {
   library.dynam.unload('rdirdirgamma', libpath)
}
