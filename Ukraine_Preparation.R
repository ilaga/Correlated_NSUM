
library(readxl)

Ukraine = read_excel("./Data/ScaleUp2009Eng.xlsx", 
                      sheet = "data")
ukraine= as.data.frame(Ukraine)
Ukraine = ukraine

# Delete unnecessary columns
# Focus on unknown groups IDU, FSW, MSW, and MSM
Ukraine = Ukraine[,c("V6", "V8", "V10", "V12", "V14", "V16", "V18", "V20", "V22",
                     "V24", "V26", "V28", "V30", "V32", "V34", "V36", "V38",
                     "V40", "V84", "V54", "V68", "V74", "V56", "V62", "V76",
                     "V42", "V92", "V93", "V94", "V95", "V96", "V97", "V98",
                     "V99", "V100", "V101", "V102", "V103", "V104", "V105", "V106",
                     "V107", "V108", "V109", "V110", "V111", "V112", "V113", "V114",
                     "V117", "V120", "V121", "V122", "V123", "V124", "V125", "V126",
                     "V127", "V128", "V129")]







# delete observations with NAs
Ukraine[Ukraine == -99] = NA
Ukraine$V128[Ukraine$V128 == 11] = NA ## DK/NA
Ukraine$V129[Ukraine$V129 == 8] = NA ## DK/NA
Ukraine$V129[Ukraine$V129 == 9] = NA ## DK/NA
Ukraine = Ukraine[complete.cases(Ukraine), ]


##################
# Demographics

Demographics = Ukraine[,55:60]
Ukraine = Ukraine[,-c(55:60)]
# code gender as 0-1
Demographics$V124[Demographics$V124 == 2] = 0

# form 4 categories for education (basic, secondary, vocational, higher)
Demographics$EducationSecondary = 0
Demographics$EducationVocational = 0
Demographics$EducationAcademic = 0
Demographics$EducationSecondary[Demographics$V126 == 4] = 1
Demographics$EducationSecondary[Demographics$V126 == 5] = 1
Demographics$EducationSecondary[Demographics$V126 == 7] = 1
Demographics$EducationSecondary[Demographics$V126 == 8] = 1
Demographics$EducationVocational[Demographics$V126 == 6] = 1
Demographics$EducationAcademic[Demographics$V126 == 9] = 1

# code nationality: 1 Ukrainian, 0 other
Demographics$V127[Demographics$V127!=1] = 0

# employment status: 1 employed, 0 unemployed
Demographics$V128[Demographics$V128==2] = 1
Demographics$V128[Demographics$V128==3] = 1
Demographics$V128[Demographics$V128==4] = 1
Demographics$V128[Demographics$V128==5] = 1
Demographics$V128[Demographics$V128==6] = 0
Demographics$V128[Demographics$V128==7] = 0
Demographics$V128[Demographics$V128==8] = 0
Demographics$V128[Demographics$V128==9] = 0
Demographics$V128[Demographics$V128==10] = 0

# internet: code 1=1, 2=0
Demographics$V129[Demographics$V129==2] = 0


Demographics= as.matrix(Demographics)
Demographics = cbind(1, Demographics)
colnames(Demographics)[1:7] = c("Intercept", "Gender", "Age", "Education", "Nationality", "Profession", "Internet")
# delete Education column
Demographics = Demographics[,-4]

###################
## Level of Respect

LoR = Ukraine[,27:54]
Ukraine = Ukraine[,-c(27:54)]
## Remove live abroad, Internet, and HIV
LoR = LoR[,-c(26, 27, 28)]
colnames(LoR) = c("Moldovans", "Romanians", "Poles", "Jews","Roma","M20-30", "M15-17", "M70+", "F20-30", "F15-17", "F70+", "Kids", 
                  "Invalids", "Pavlo", "DivorcedMen2007", "Birth2007","Prisoners2007","IDU","MedicalDoctors","FSW", "MSW", "PhD","MSM",
                  "Nurses","Policemen")
## Add Died2007
LoR = cbind(LoR, 3)
colnames(LoR)[26] = "Died2007"
## Match ordering to Ukraine
colnames(Ukraine)= c("M20-30", "M15-17", "M70+", "F20-30", "F15-17", "F70+", "Kids", "Moldovans", "Romanians", "Poles", "Jews","Roma",
                      "Invalids", "MedicalDoctors", "Died2007", "Pavlo", "Prisoners2007", "DivorcedMen2007","Policemen", "Birth2007", "PhD",
                      "Nurses","FSW", "MSW", "MSM", "IDU")
LoR = LoR[,match(colnames(Ukraine), colnames(LoR))]
LoR[LoR==8] = 3 ## Set missing LoR to 3
mLoR= colMeans(LoR)

# Upper truncate values at 150 (some values are like 15000)
Ukraine[Ukraine > 150] = 150


# Convert data to a matrix
Ukraine = as.matrix(Ukraine)

# Save the data

write.table(Ukraine, "./Data/Ukraine_y.csv", sep=",", row.names = F, col.names = F)
write.table(Demographics, "./Data/Ukraine_x_cov.csv", sep=",", row.names = F, col.names = F)
write.table(LoR, "./Data/resplevel_Ukraine_z_cov.csv", sep=",", row.names = F, col.names = F)





# Create plots for Figure 1 -----------------------------------------------

library(ggplot2)

y = as.matrix(read.csv("./Data/Ukraine_y.csv", header = F))
known <- c(4088438, 935153, 1328606, 3966956, 889443, 2993846, 1869098, 258619, 150989, 144130,104186, 47587, 278195, 222884, 762877,
           323863, 108511, 178364, 273200, 472657, 69471, 487148)
names(known)<- c("M20-30", "M15-17", "M70+", "F20-30", "F15-17", "F70+", "Kids", "Moldovans", "Romanians", "Poles", "Jews","Roma",
                 "Invalids", "MedicalDoctors", "Died2007", "Pavlo", "Prisoners2007", "DivorcedMen2007","Policemen", "Birth2007", "PhD", "Nurses")


colnames(y) = c(names(known), "FSW", "MSW", "MSM", "IDU")
y[y > 30] = 30 ## Truncate for visualization only
gg.y = reshape2::melt(y)


gg.f70 = ggplot(subset(gg.y, Var2 == "F70+")) + geom_bar(aes(x = value)) +
  ggtitle("F70+") + xlab("Response") + ylab("Frequency") +
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.x = element_text(size = 25, face = "bold"),
        axis.title.y = element_text(size = 25, face = "bold"),
        title = element_text(size = 25, face = "bold"))

gg2.kids = ggplot(subset(gg.y, Var2 == "Kids")) + geom_bar(aes(x = value)) +
  ggtitle("Kids") + xlab("Response") + ylab("Frequency") +
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.x = element_text(size = 25, face = "bold"),
        axis.title.y = element_text(size = 25, face = "bold"),
        title = element_text(size = 25, face = "bold"))

gg3.idu = ggplot(subset(gg.y, Var2 == "IDU")) + geom_bar(aes(x = value)) +
  ggtitle("IDU") + xlab("Response") + ylab("Frequency") +
  theme(axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.x = element_text(size = 25, face = "bold"),
        axis.title.y = element_text(size = 25, face = "bold"),
        title = element_text(size = 25, face = "bold"))



ggsave("./Figures/Ukraine_F70_barplot.jpg", gg.f70,
       width = 6, height = 6)

ggsave("./Figures/Ukraine_Kids_barplot.jpg", gg2.kids,
       width = 6, height = 6)

ggsave("./Figures/Ukraine_IDU_barplot.jpg", gg3.idu,
       width = 6, height = 6)
