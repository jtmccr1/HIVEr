accuracy_stringent<-read.csv("../../Analysis/Host_level_IAV_evolution/data/reference/accuracy_stringent.csv",
                             stringsAsFactors = F)
accuracy_stringent<-dplyr::select(accuracy_stringent,-dplyr::starts_with("X"))

devtools::use_data(accuracy_stringent,overwrite = T)
