##### Relationship between carriage duration and invasiveness ####

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)

#### data prep ####
carrdur <- read.csv("Q:/Technical/Python/seroparams.csv")
carrdur <- carrdur %>% select(serotype, maxcarrdur)
colnames(carrdur) <- c('Serotype', 'maxcarrdur')
carrdur <- carrdur[which(carrdur$Serotype != '15B/C_old'),]
carrdur <- carrdur[which(carrdur$Serotype != '15B'),]
carrdur <- carrdur[which(carrdur$Serotype != '15C'),]
VT7_index <- which(carrdur$Serotype %in% c("4", "6B", "9V", "14", "18C", "19F", "23F"))
VT10_index <- which(carrdur$Serotype %in% c("1", "5", "7F"))
VT13_index <- which(carrdur$Serotype %in% c("3", "6A", "19A"))
VT_index <- c(VT7_index, VT10_index, VT13_index)
carrdur$Serogroup <- rep("NVT", nrow(carrdur))
carrdur$Serogroup[VT7_index] <- "VT7"
carrdur$Serogroup[VT10_index] <- "VT10"
carrdur$Serogroup[VT13_index] <- "VT13"

fi_inv <- read.csv("Q:/Technical/Python/finlandlocalinv.csv")
colnames(fi_inv) <- c('Serotype', 'localinv', 'localinv_low', 'localinv_high')
VT7_index <- which(fi_inv$Serotype %in% c("4", "6B", "9V", "14", "18C", "19F", "23F"))
VT10_index <- which(fi_inv$Serotype %in% c("1", "5", "7F"))
VT13_index <- which(fi_inv$Serotype %in% c("3", "6A", "19A"))
VT_index <- c(VT7_index, VT10_index, VT13_index)
fi_inv$Serogroup <- rep("NVT", nrow(fi_inv))
fi_inv$Serogroup[VT7_index] <- "VT7"
fi_inv$Serogroup[VT10_index] <- "VT10"
fi_inv$Serogroup[VT13_index] <- "VT13"

fi_inv <- read.csv("Q:/Technical/Python/francelocalinv.csv")
colnames(fi_inv) <- c('Serotype', 'localinv', 'localinv_low', 'localinv_high')
VT7_index <- which(fi_inv$Serotype %in% c("4", "6B", "9V", "14", "18C", "19F", "23F"))
VT10_index <- which(fi_inv$Serotype %in% c("1", "5", "7F"))
VT13_index <- which(fi_inv$Serotype %in% c("3", "6A", "19A"))
VT_index <- c(VT7_index, VT10_index, VT13_index)
fi_inv$Serogroup <- rep("NVT", nrow(fi_inv))
fi_inv$Serogroup[VT7_index] <- "VT7"
fi_inv$Serogroup[VT10_index] <- "VT10"
fi_inv$Serogroup[VT13_index] <- "VT13"

globalinv <- read.csv("Q:/Technical/R/Case-to-carrier/consol.child-new.csv")
globalinv$X <- NULL
colnames(globalinv) <- c('Serotype', 'globalinv', 'globalinv_low', 'globalinv_high', 'Serogroup')

localdf <- left_join(carrdur, fi_inv)
alldf <- left_join(localdf, globalinv)

cols <- c("VT7" = "#35B779FF", "VT10"= "#31688EFF", "VT13" = "#E69F00", "NVT" = "#440154FF")

ggplot(alldf) +
  geom_point(aes(x = maxcarrdur, y = localinv, group = Serogroup, colour = Serogroup)) +
  geom_errorbar(aes(x = maxcarrdur, ymin = localinv_low, ymax = localinv_high, colour = Serogroup)) +
  ggrepel::geom_text_repel(aes(x = maxcarrdur, y = localinv, label = Serotype, group = Serogroup, 
                             colour = Serogroup), fontface = "bold",show.legend = FALSE) +
  scale_color_manual(values= cols) + scale_y_continuous(trans = 'log10') +
  labs(x = 'Max carriage duration (days)', y = 'Finland local invasiveness', colour = 'Category')
