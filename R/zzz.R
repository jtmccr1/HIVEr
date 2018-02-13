.onAttach <- function(pkgname, libname) {
  my_palette <- wesanderson::wes_palette("Zissou")
  my_theme <- ggplot2::theme_classic()+ ggplot2::theme(text=ggplot2::element_text(family="Arial",size = 18))
  ggplot2::theme_set(my_theme)

  assign("scale_colour_discrete", function(..., values = my_palette) scale_colour_manual(..., values = values), globalenv())
  assign("scale_fill_discrete", function(..., values = my_palette) scale_fill_manual(..., values = values), globalenv())
}
