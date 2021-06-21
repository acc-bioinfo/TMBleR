# data(ExamplePaneldesign)
# data(vcfs_NoCancer_ForPanel)
# data(HumanTumorSuppressors)
# data(Hellman_SimulatedFM1Panel_WES)
# data(COSMIC_hg19)
# data(COSMIC_hg38)
utils::globalVariables(c( "vcfs_NoCancer_ForPanel"
                        , "HumanTumorSuppressors"
                        , "Hellman_SimulatedFM1Panel_WES"
                        , "COSMIC_hg19_demo"
                        , "COSMIC_hg38_demo"
                        , "COSMIC_hg19"
                        , "COSMIC_hg38"
                        , ".")
                       , package = "TMBleR", add = TRUE)

# Set a debugging variable
debug_env <- new.env()
#debug_env$debug = FALSE
debug_env$debug = FALSE