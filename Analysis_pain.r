######################################################################################
#####  R script for item content overlap in central sensitization questionnaires  ####                                                                         ##### 
#####  Based on Fried (2016;2018)                                                 ####
#####                                                                             #### 
######################################################################################

setwd("M:/Documents/2023_Content_Pain_Sensitization/R_code_pain")


library('qgraph')
library('ggplot2')
library('data.table')
library('reshape2')
library('psych')
library('ade4')
library('viridis')
library('geomtextpath')
library('dplyr')


##### Data preparation
data = fread("MatrixB.csv")     #Load data for estimating Jaccard index (no difference between specific and compound symptoms)
dataplot = fread("MatrixA.csv") #Load data for plot (difference between specific and compound symptoms)


##### Estimation of overlap, using the Jaccard Index

# GPQ
data1<-data[which(data$painDETECT==1|data$GPQ==1),]
a1<-1-(dist.binary(matrix(c(data$GPQ ,data$PSQ),nrow=2,byrow=T), method = 1)^2); a1
b1<-1-(dist.binary(matrix(c(data$GPQ ,data$CSI),nrow=2,byrow=T), method = 1)^2)
c1<-1-(dist.binary(matrix(c(data$GPQ ,data$CAP_Knee),nrow=2,byrow=T), method = 1)^2)
d1<-1-(dist.binary(matrix(c(data$GPQ ,data$painDETECT),nrow=2,byrow=T), method = 1)^2)
GPQ.v<-c(1,a1,b1,c1,d1)

# PSQ
a2<-1-(dist.binary(matrix(c(data$PSQ, data$GPQ),nrow=2,byrow=T), method = 1)^2)
b2<-1-(dist.binary(matrix(c(data$PSQ ,data$CSI),nrow=2,byrow=T), method = 1)^2)
c2<-1-(dist.binary(matrix(c(data$PSQ ,data$CAP_Knee),nrow=2,byrow=T), method = 1)^2)
d2<-1-(dist.binary(matrix(c(data$PSQ ,data$painDETECT),nrow=2,byrow=T), method = 1)^2)
PSQ.v<-c(a2,1,b2,c2,d2)

# CSI
a3<-1-(dist.binary(matrix(c(data$CSI, data$GPQ),nrow=2,byrow=T), method = 1)^2)
b3<-1-(dist.binary(matrix(c(data$CSI, data$PSQ),nrow=2,byrow=T), method = 1)^2)
c3<-1-(dist.binary(matrix(c(data$CSI, data$CAP_Knee),nrow=2,byrow=T), method = 1)^2)
d3<-1-(dist.binary(matrix(c(data$CSI, data$painDETECT),nrow=2,byrow=T), method = 1)^2)
CSI.v<-c(a3,b3,1,c3,d3)

# CAP-Knee
a4<-1-(dist.binary(matrix(c(data$CAP_Knee, data$GPQ),nrow=2,byrow=T), method = 1)^2)
b4<-1-(dist.binary(matrix(c(data$CAP_Knee, data$PSQ),nrow=2,byrow=T), method = 1)^2)
c4<-1-(dist.binary(matrix(c(data$CAP_Knee, data$CSI),nrow=2,byrow=T), method = 1)^2)
d4<-1-(dist.binary(matrix(c(data$CAP_Knee, data$painDETECT),nrow=2,byrow=T), method = 1)^2)
CAP_Knee.v<-c(a4,b4,c4,1,d4)

# painDETECT
a5<-1-(dist.binary(matrix(c(data$painDETECT, data$GPQ),nrow=2,byrow=T), method = 1)^2)
b5<-1-(dist.binary(matrix(c(data$painDETECT, data$PSQ),nrow=2,byrow=T), method = 1)^2)
c5<-1-(dist.binary(matrix(c(data$painDETECT, data$CSI),nrow=2,byrow=T), method = 1)^2)
d5<-1-(dist.binary(matrix(c(data$painDETECT, data$CAP_Knee),nrow=2,byrow=T), method = 1)^2)
painDETECT.v<-c(a5,b5,c5, d5,1)



# Create table
M = matrix(nrow=5, ncol=5) 
colnames(M) <- c("GPQ",	"PSQ",	"CSI",	"CAP-Knee",	"painDETECT")
rownames(M) <- c("GPQ",	"PSQ",	"CSI",	"CAP-Knee",	"painDETECT")
M[1,]<-GPQ.v
M[2,]<-PSQ.v
M[3,]<-CSI.v
M[4,]<-CAP_Knee.v
M[5,]<-painDETECT.v
isSymmetric(M)
M



M[M == 1] <- 0 # replace diagonal with 0
colMeans(M)
mean(colMeans(M)) # 0.11, is underestimation, see Fried (2018)



# Correct mean overlap 
x <- colMeans(M)
x[1] <- sum(M[,1])/4
x[2] <- sum(M[,2])/4
x[3] <- sum(M[,3])/4
x[4] <- sum(M[,4])/4
x[5] <- sum(M[,5])/4

mean(x) # 0.14 without the diagonal, which is correct


length1<-c(7,17,25,8,10) # length of original questionnaires; e.g. PSQ 17 items originally
length2<-c(7,3,25,8,10) # items in analysis per scale; PSQ captures 3 items

cor(length1, colMeans(M)) # -0.62
cor(length2, colMeans(M)) # -0.84



##### Compute mean symptom appearance across questionnaires
# Sum across rows to get the number of "1" occurrences per row
row_sums = rowSums(data)

# Calculate the mean number of "1" occurrences
mean_ones = mean(row_sums)
# Calculate the median number of "1" occurrences
median_ones = median(row_sums)
# Calculate the mode number of "1" occurrences
mode_ones = mode(row_sums)

# Function to calculate mode
getMode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Sum across rows to get the number of occurrences per row
row_sums = rowSums(data)
# Calculate the mode of the number of occurrences
mode_ones = getMode(row_sums)

# Print the results
print(paste("Mean number of appearances: ", mean_ones))
print(paste("Median number of appearances: ", median_ones))
print(paste("Mode number of appearances: ", mode_ones))


# Generate frequency table
freq_table <- table(row_sums)
# Calculate percentages
percentages <- prop.table(freq_table) * 100
# Combine frequencies and percentages in a data frame for better readability
freq_table_df <- data.frame(Number_of_1s = names(freq_table), Frequency = freq_table, Percentage = percentages)

# Print the frequency table with percentages
print(freq_table_df)



##### Figure 1 
set.seed(223)
d <- dataplot
d[, S := factor(paste0("S",1:nrow(d)))] #Create symptom variable
d = data.table::melt(d, id.vars="S", variable.name="Scales", value.name="Type") #Transform to long format
d = d[Type>=1] #Keep the scales in which the symptoms are 1 (present) or 2 (included)
d[, Type := factor(Type, labels=c("Hypothetical symptom", "Actual symptom"))]
d[, count := .N, by=S]


# Add a column with with the higher order categories.
d <- d %>% 
  mutate(category = case_when(
    S %in% c('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7') ~ 'Nociplastic pain manifestations',
    S %in% c('S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19',
             'S20', 'S21', 'S22', 'S23', 'S24', 'S25', 'S26', 'S27', 'S28', 'S29', 'S30',
             'S31', 'S32', 'S33') ~ 'Psychosomatic symptoms',
    S == 'S34' ~ 'Pain location',
    S == 'S35' ~ 'Pain course',
    S %in% c('S36', 'S37', 'S38', 'S39') ~ 'Neuropathic pain descriptors',
    TRUE ~ NA_character_ # for any other case not covered above
  ))


# Symptom order
sympt.order = d[, .N, by=S][order(N)][, S] #Replace by order
d[, S := factor(S, levels = sympt.order)]

# Scale order by frequency
scale.order = d[, .N, by=Scales][order(N)][, Scales]
d[, Scales := factor(Scales, levels = scale.order)]
d[, Scales2 := as.numeric(Scales)]




# Reordering the d dataframe for the flower plot

# First, convert 'S' and 'Scales' to character to avoid factor issues
d[, S := as.character(S)]
d[, Scales := as.character(Scales)]
d[, category := as.character(category)] # Make sure category is also character

# Create a combined key of 'category' and 'S' to ensure uniqueness
d[, combined_key := paste(category, S, sep = "_")]

# Order the combined keys by their frequency
combined_order <- d[, .N, by = combined_key][order(-N), combined_key]

# Assign the ordered combined keys as factor levels to 'S'
d[, S := factor(combined_key, levels = combined_order)]

# Separate the 'S' and 'category' again
d[, c("category", "S") := tstrsplit(S, "_", fixed = TRUE)]
# Convert 'S' back to factor with the levels in the order of appearance
d[, S := factor(S, levels = unique(S))]

# Do the same for 'Scales'
d[, combined_key := paste(category, Scales, sep = "_")]
combined_scale_order <- d[, .N, by = combined_key][order(-N), combined_key]
d[, Scales := factor(combined_key, levels = combined_scale_order)]
d[, Scales := tstrsplit(Scales, "_", fixed = TRUE)[[2]]]
d[, Scales := factor(Scales, levels = unique(Scales))]

# Create a numeric representation of 'Scales'
d[, Scales2 := as.numeric(Scales)]

# Now remove the temporary 'combined_key'
d[, combined_key := NULL]

# Check the result
View(d)









# Plot
pal1 <- c("#D55E00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")



# Create a named vector for new labels
new_labels <- c("GPQ" = "GPQ",
                "PSQ" = "PSQ",
                "CSI" = "CSI", 
                "CAP" = "CAP-Knee", 
                "painDETECT" = "painDETECT")


plot_flower <- ggplot(d, aes(x=S, y=Scales2, group=S, color=as.factor(Scales), 
                             shape=Type, rev=F)) +
  geom_line() + # keep this here, otherwise there is an error
  xlab("") +
  ylab("") +
  geom_hline(yintercept = 1:5, colour = "grey80", size = .2) +
  geom_vline(xintercept = 1:39, colour = "grey80", size = .2) +
  geom_line(colour="grey60") +
  geom_point(size=2.8, fill="white", stroke=1.5) + # Increase stroke width
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=.6, fill="white", color=NA) +
  coord_curvedpolar() +
  scale_shape_manual(values=c(21,19)) +
  theme(axis.text.x = element_text(size = 9.5, vjust = 0.5),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "left",
        legend.title = element_blank(),
        plot.margin = unit(rep(.5,4), "lines")) +
  labs(fill="") +
  scale_y_continuous(limits=c(-3,5), expand=c(0,0), breaks=1:8, labels=NULL) +
  scale_color_manual(values=pal1) + 
  scale_color_manual(values=pal1, labels=new_labels) +
  guides(shape = guide_legend(position = "bottom"))

plot_flower

ggsave(plot=plot_flower,filename="Figure1.pdf", width=7, height=6, useDingbats=FALSE)
# Figure 1 was further adjusted in Inkscape (e.g., item legend)



ggsave("plot.tiff", width = 2187, height = 1351, units = "px", dpi=300)


