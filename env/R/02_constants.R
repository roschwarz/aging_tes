volcano_colors <- setNames(
    c("#e63946", "gray", "#457b9d"),
    c("Up", "Not Sig", "Down")
)

chr_of_interest <- paste0("chr", c(1:19, "X", "Y"))

FDR <- 0.05

tissues = c("brain", "skin", "blood")
orders.of.interest <- c("LINE", "SINE", "LTR", "DNA")


# set the theme of interest
#theme_set(theme_rob(12))

# color codes for specific scenarios
#tissue.color <- c(brain = '#264653', skin = '#2A9D8F',  blood = '#E9C46A')
tissue.color <- c(background = 'black',brain = '#58B2AA', skin = '#EFA081',  blood = '#685299')
tissue.text.color <- c( brain = "#ffffff", skin = "#000000", blood = "#ffffff")

#order.color = c(DNA = "#798E87", LINE = "#C27D38", LTR = "#CCC591", SINE = "#29211F")
order.color = c(DNA = "#E49400", LINE = "#831335", LTR = "#14273F", SINE = "#8493AE")
#direction.color = c(up = '#d1495b', down = '#66a182', equal = 'grey' )
direction.color = c(up = '#E63846', down = '#467B9D', equal = 'grey80' )