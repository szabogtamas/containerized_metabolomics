setHook("rstudio.sessionInit", function(newSession) {
  if (newSession)
    message("Hi, this is a fork of Rocker project.")
    message("Welcome to RStudio ", rstudioapi::getVersion())
    rstudioapi::applyTheme("Solarized Dark")
}, action = "append")