library(shiny)
library(ggnewscale)
library(ggplot2)
library(glmmTMB)
library(nnet)

source('/Users/sophie/Documents/implied_weights/impliedweights_randeffs.R')

## Import data
#study = read.csv('/Users/sophie/Documents/SpatialConf/archived/Study_dataset_2010.csv')
load('../mapdat.RData') # created by preprocessing.R

#study$FIPS = str_pad(study$FIPS, 5, pad = '0')
#study$pm25_binary = ifelse(study$qd_mean_pm25 > 9, 1, 0) # naaqs standard
# df$pm25_binary = ifelse(df$mean_pm25 > 9, 1, 0) # naaqs standard
# qs = quantile(df$mean_pm25, probs = c(0.2, 0.4, 0.6, 0.8))
# qs = c(-Inf, qs, Inf)
# # Bin PM2.5 by qs
# df$pm25_multi = cut(df$mean_pm25, breaks = qs, labels = FALSE)
# 
# # Create design matrix. 
# df$statefactor = as.numeric(as.factor(df$STATEFP10))
#studyexclude = which(study$STATE_CODE %in% c("CA", "OR", "WA", "ID", "NV", "UT", "AZ", "CO", "NM", "WY", "MT"))
#study = study[-studyexclude,]

# study = mutate(study,
#                STATEFP10 = str_sub(FIPS, 1, 2),
#                COUNTYFP10 = str_sub(FIPS, 3, 5))

# Import adjacency matrix created from spacebench script
#adj = read.csv("/Users/sophie/Documents/SpatialConf/archived/adjacency_matrix.csv", header = F)

# column bind adj to study
#study = cbind(study, adj[-studyexclude,-studyexclude])

#states = st_read("/Users/sophie/Documents/SpatialConf/archived/tl_2010_us_county10/tl_2010_us_county10.shp")
stateboundaries = st_read("/Users/sophie/Downloads/tl_2010_us_state10/tl_2010_us_state10.shp")
# 
# # Omit hispanic (collinear with other races)
# #mapdat = right_join(states, study, by = c("STATEFP10", "COUNTYFP10"))
# # Order mapdat by statefactor
# mapdat = df
# mapdat = mapdat[order(mapdat$statefactor),]
clusters = mapdat$statefactor
nks = table(clusters)

# Extract A as the last 2695 columns of mapdat
#A = st_drop_geometry(mapdat[,(ncol(mapdat)-2695):(ncol(mapdat)-1)])
#A = as.matrix(A)
X = st_drop_geometry(mapdat[,20:32])
X = cbind(1, X)
X = as.matrix(X)
Z = mapdat$pm25_binary
Y = mapdat$all_all_65up_imputed
state.dummy = model.matrix(~STATEFP10 - 1, data = mapdat)
colnames(state.dummy)=gsub("STATEFP10","",as.character(colnames(state.dummy)))
Xstate = cbind(1, st_drop_geometry(mapdat[,20:32]), state.dummy[,-1]) 
Xstate = as.matrix(Xstate)

# Compute fixed effects weights
mapdat$wts_fe = fixedeffects(
  X = Xstate,
  Z = Z
)
# Compute complete pooling weights
mapdat$wts_pooled = fixedeffects(
  X = X,
  Z = Z
)

# find states whose counties are only all 0 or all 1
for (st in unique(mapdat$STATEFP10)){
  if (length(unique(mapdat$pm25_binary[mapdat$STATEFP10 == st])) == 1){
    print(c(st, sum(mapdat$wts_fe[mapdat$STATEFP10 == st])))
  }
}

# Load in CAR and SAR weights (precomputed)
load("../carweights702.rdata")
load("../sarweights702.rdata")
mapdat$carweights = carweights
mapdat$sarweights = sarweights
load("../carweights702_binary.rdata")
load("../sarweights702_binary.rdata")
mapdat$carweights_binary = carweights
summary(mapdat$carweights_binary)
mapdat$sarweights_binary = sarweights

p = ncol(X)
K = length(unique(clusters))
# Xks = lapply(1:K, function(k) X[clusters == k,])
# Zks = lapply(1:K, function(k) Z[clusters == k])
# nks = sapply(1:K, function(k) sum(clusters == k))
# mks = sapply(1:K, function(k) sum(Z[clusters == k]))

# Xbarks = lapply(1:K, function(k) {
#   Xk = matrix(X[clusters == k,], ncol = p, nrow = nks[k])
#   return(matrix(colMeans(Xk), ncol = p, nrow = 1))
# }
# )

# Xtbarks = lapply(1:K, function(k) {
#   if (mks[k] > 0){
#     if (nks[k] == 1){
#       return(Xks[[k]])
#     }
#     else{
#       Xktreated = matrix(Xks[[k]][Zks[[k]] == 1,], nrow = mks[k], ncol = p)
#       return(matrix(colMeans(Xktreated), ncol = p, nrow = 1))
#     }
#   }
#   else{
#     return(matrix(0, ncol = p, nrow = 1))
#   }
# }
# )


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel(h1("The Trinity Explained: Implied weights of Pooled, Random Effect and Fixed Effect models", align = "center")),
  
  sidebarPanel(
    # Input: Slider for sig2gam ----
    sliderInput(inputId = "sig2gam",
                label = "sig2gam:",
                min = 0,
                max = 1,
                value = 0.6),
    # Input: Slider for sig2eps ----
    sliderInput(inputId = "sig2eps",
                label = "sig2eps:",
                min = 0.00001,
                max = 10,
                value = 7)#,
    # Input: slider for phi ----
    # sliderInput(inputId = "phi",
    #             label = "phi:",
    #             min = 0,
    #             max = 0.999,
    #             value = 0.5),
    #checkboxInput("fe_weights", "Show Fixed Effects implied weights", TRUE)
  ),
  
  tabsetPanel(
    tabPanel("pooled/RE/FE",
             fluidRow(
               column(4, #style='padding:5px;',
                      h4("Complete Pooling", align = "center"),
                      plotOutput(outputId = "pooledPlot")),
               column(4, #style='padding:5px;',      
                      h4("Random Effects", align = "center"),
                      plotOutput(outputId = "rePlot")),
               column(4, #style='padding:5px;',
                      h4("Fixed Effects", align = "center"),
                      plotOutput(outputId = "fePlot"))
             )
             ),
    tabPanel("Treated/Control",
             fluidRow(
               column(6,
                      h4("Binary Treatment", align = "center"),
                      plotOutput(outputId = "binaryPlot")),
               column(6,
                      h4("Multivalued Treatment", align = "center"),
                      plotOutput(outputId = "multiPlot"))
             )
             ),
    tabPanel("CAR/SAR",
             fluidRow(
               column(6,
                      h4("CAR Weights", align = "center"),
                      plotOutput(outputId = "carPlot")),
               column(6,
                      h4("SAR Weights", align = "center"),
                      plotOutput(outputId = "sarPlot"))
             )
    ),
    tabPanel("CAR/SAR Binary Adjmat",
             fluidRow(
               column(6,
                      h4("CAR Weights", align = "center"),
                      plotOutput(outputId = "carPlot_binary")),
               column(6,
                      h4("SAR Weights", align = "center"),
                      plotOutput(outputId = "sarPlot_binary"))
             )
    ),
    tabPanel("Balance",
             fluidRow(
               column(6,
                      align = "center",
                      h4("Covariate ASMD", align = "center"),
                      tableOutput("balance")),
               column(6,
                      align = "center",
                      h4("State-Indicator ASMD", align = "center"),
                      tableOutput("statebalance"))#,
               # column(2,
               #        align = "center",
               #        h4("Empirical Bayes (x 10^4)", align = "center"),
               #        tableOutput("empiricalBayes"))
             )
            ),
    tabPanel("Dispersion",
             fluidRow(
               column(4,
                      align = "center",
                      h4("Dispersion of Weights", align = "center"),
                      tableOutput("dispersion")
                      ),
               column(8,
                      align = "center",
                      h4("Histograms of Implied Weights", align = "center"),
                      plotOutput(outputId = "weighthists")
               )
              )
    ),
    tabPanel("Outcome",
             fluidRow(
               column(6,
                      h4("Heart Disease Hospitalization Rate 2016-2018", align = "center"),
                      plotOutput(outputId = "outcomeplot")),
               column(6,
                      h4("Heart Disease Hospitalization Rate 2016-2018, Imputed", align = "center"),
                      plotOutput(outputId = "outcome_imputed_plot"))
             )
    ),
    tabPanel("Estimates",
             fluidRow(
               column(12,
                      align = "center",
                              h4("Treatment Effect Estimates", align = "center"),
                              tableOutput("estimates"))
               ),
             fluidRow(
               column(12,
                      align = "center",
                      h4("Multitreatment Effect Estimates", align = "center"),
                      tableOutput("multiestimates"))
             )
             ),
    tabPanel("Multivalued Treatment",
             img(src='IW_treatment_Z2.png', height="75%", width="75%", align = "center"),
             img(src='IW_treatment_Z3.png', height="75%", width="75%", align = "center"),
             img(src='IW_treatment_Z4.png', height="75%", width="75%", align = "center"),
             img(src='IW_treatment_Z5.png', height="75%", width="75%", align = "center"))
             #,
    # tabPanel("Within-State Balance (RE only)",
    #          tableOutput("balancebystate"),
    #          tableOutput("balancebystate_QX"))
  )
  
    
    # Sidebar panel for inputs ----
    # fluidRow(
    #   column(6,
    # 
    #   # Input: Slider for sig2gam ----
    #   sliderInput(inputId = "sig2gam",
    #               label = "sig2gam:",
    #               min = 0,
    #               max = 1,
    #               value = 0.5),
    #   # Input: Slider for sig2eps ----
    #   sliderInput(inputId = "sig2eps",
    #               label = "sig2eps:",
    #               min = 0.00001,
    #               max = 10,
    #               value = 7),
    #   checkboxInput("fe_weights", "Show Fixed Effects implied weights", TRUE)
    #   ),
      # column(6,
      #        tableOutput("balance")
      #        )#,
      # column(4,
      #        textOutput("todo")
      #        )
    
    
    # Main panel for displaying outputs ----
    # mainPanel(
    #   plotOutput(outputId = "rePlot"),
    #   plotOutput(outputId = "fePlot")
    # )
    
    
    # Create fluidRow below for sliders and checkbox
    # fluidRow(
    #   column(6, #style='padding:5px;',
    #          wellPanel(
    #            sliderInput(inputId = "sig2gam",
    #                        label = "sig2gam:",
    #                        min = 0,
    #                        max = 10,
    #                        value = 0.5),
    #            br(),
    #            sliderInput(inputId = "sig2eps",
    #                        label = "sig2eps:",
    #                        min = 0.00001,
    #                        max = 10,
    #                        value = 7)
    #            )
    #          ),
    #   column(6, #style='padding:5px;',
    #          checkboxInput("fe_weights", "Show Fixed Effects implied weights", TRUE))
    # )
  )

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  global = reactiveValues(balance = NULL)
  
  wts_re = reactive(as.vector(
    impliedweightsslow(
      X = X,
      Z = Z,
      sig2gam = input$sig2gam,
      sig2eps = input$sig2eps,
      clusters = clusters
    )
  ))
  
  # carweights = reactive(as.vector(
  #   car(X = X,
  #       Z = Z,
  #       phi = input$phi,
  #       A = A,
  #       sig2eps = sig2eps
  #   )
  # ))
  # 
  # sarweights = reactive(as.vector(
  #   sar(X = X,
  #       Z = Z,
  #       phi = input$phi,
  #       A = A,
  #       sig2eps = sig2eps
  #   )
  # ))
  
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs change
  # 2. Its output type is a plot
  output$rePlot <- renderPlot({
    
    # wts_re = impliedweightsfast(X = X,
    #                             Xks = Xks,
    #                             Zks = Zks,
    #                             nks = nks,
    #                             mks = mks,
    #                             Xbarks = Xbarks,
    #                             Xtbarks = Xtbarks,
    #                             sig2gam = input$sig2gam,
    #                             sig2eps = input$sig2eps)
    mapdat$wts_re = wts_re()
    ggplot(mapdat) +
      xlim(-125, -65) +
      ylim(25, 50) +
      geom_sf(aes(fill = wts_re), color=NA, size = 0.005) +
      geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
      theme_minimal() +
      scale_fill_gradient2(expression('weights'),
                           low = "#8b0000",
                           mid = "white",
                           high = "navy",
                           midpoint = 0,
                           breaks = c(-0.001, -0.0005, 0.0005, 0.001, 0.002),
                           labels = c("-0.001", "-0.0005", "0.0005", "0.001", "0.002"),
                           na.value = "grey",
                           limits = c(-0.002, 0.004)) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            line = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text.align = 0.75,
            legend.key.width = unit(100, "points"),
            panel.grid.major = element_line(colour = "transparent"))
    # ggplot() +
    #   xlim(-125, -65) +
    #   ylim(25, 50) +
    #   geom_sf(data = mapdat[mapdat$pm25_binary == 0, ], aes(fill = wts_re), color=NA, size = 0.005) +
    #   scale_fill_gradient2(expression('weights'),
    #                       low = "white",
    #                       mid = "white",
    #                       high = "#8b0000",
    #                       midpoint = 0,
    #                       breaks = c(-0.001, -0.0005, 0.002, 0.003, 0.004),
    #                       labels = c("-0.001", "-0.0005", "0.002", "0.003", "0.004"),
    #                       na.value = "grey",
    #                       limits = c(-0.002, 0.004)) +
    #   new_scale_fill() +
    #   geom_sf(data = mapdat[mapdat$pm25_binary == 1, ], aes(fill = wts_re), color=NA, size = 0.005) +
    #   scale_fill_gradient2(expression('weights'),
    #                        low = "white",
    #                        mid = "white",
    #                        high = "navy",
    #                        midpoint = 0,
    #                        breaks = c(-0.001, -0.0005, 0.002, 0.003, 0.004),
    #                        labels = c("-0.001", "-0.0005", "0.002", "0.003", "0.004"),
    #                        na.value = "grey",
    #                        limits = c(-0.002, 0.004)) +
    #   geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
    #   theme_minimal() +
    #   theme(
    #     axis.text.x = element_blank(),
    #     axis.text.y = element_blank(),
    #     axis.ticks = element_blank(),
    #     line = element_blank(),
    #     axis.title = element_blank(),
    #     legend.position = "bottom",
    #     legend.direction = "horizontal",
    #     legend.text.align = 0.75,
    #     legend.key.width = unit(100, "points"),
    #     panel.grid.major = element_line(colour = "transparent")
    #   ) +
    #   guides(fill = guide_colourbar(barwidth = 5, title.position = "top", label.position = "bottom"))
  })
  
  output$fePlot <- renderPlot({
        ggplot(mapdat) +
        xlim(-125, -65) +
        ylim(25, 50) +
        geom_sf(aes(fill = wts_fe), color=NA, size = 0.005) +
        geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
        theme_minimal() +
        scale_fill_gradient2(expression('weights'),
                             low = "#8b0000",
                             mid = "white",
                             high = "navy",
                             midpoint = 0,
                             breaks = c(-0.001, -0.0005, 0.0005, 0.001, 0.002),
                             labels = c("-0.001", "-0.0005", "0.0005", "0.001", "0.002"),
                             na.value = "grey",
                             limits = c(-0.002, 0.004)) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              line = element_blank(),
              axis.title = element_blank(),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text.align = 0.75,
              legend.key.width = unit(100, "points"),
              panel.grid.major = element_line(colour = "transparent"))
    # ggplot() +
    #   xlim(-125, -65) +
    #   ylim(25, 50) +
    #   geom_sf(data = mapdat[mapdat$pm25_binary == 0, ], aes(fill = wts_fe), color=NA, size = 0.005) +
    #   scale_fill_gradient2(expression('weights'),
    #                        low = "white",
    #                        mid = "white",
    #                        high = "#8b0000",
    #                        midpoint = 0,
    #                        breaks = c(-0.001, -0.0005, 0.002, 0.003, 0.004),
    #                        labels = c("-0.001", "-0.0005", "0.002", "0.003", "0.004"),
    #                        na.value = "grey",
    #                        limits = c(-0.002, 0.004)) +
    #   new_scale_fill() +
    #   geom_sf(data = mapdat[mapdat$pm25_binary == 1, ], aes(fill = wts_fe), color=NA, size = 0.005) +
    #   scale_fill_gradient2(expression('weights'),
    #                        low = "white",
    #                        mid = "white",
    #                        high = "navy",
    #                        midpoint = 0,
    #                        breaks = c(-0.001, -0.0005, 0.002, 0.003, 0.004),
    #                        labels = c("-0.001", "-0.0005", "0.002", "0.003", "0.004"),
    #                        na.value = "grey",
    #                        limits = c(-0.002, 0.004)) +
    #   geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
    #   theme_minimal() +
    #   theme(
    #     axis.text.x = element_blank(),
    #     axis.text.y = element_blank(),
    #     axis.ticks = element_blank(),
    #     line = element_blank(),
    #     axis.title = element_blank(),
    #     legend.position = "bottom",
    #     legend.direction = "horizontal",
    #     legend.text.align = 0.75,
    #     legend.key.width = unit(100, "points"),
    #     panel.grid.major = element_line(colour = "transparent")
    #   ) #+
      #guides(fill = guide_colourbar(barwidth = 1))
  })
  
  output$pooledPlot <- renderPlot({
    ggplot(mapdat) +
      xlim(-125, -65) +
      ylim(25, 50) +
      geom_sf(aes(fill = wts_pooled), color=NA, size = 0.005) +
      geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
      theme_minimal() +
      scale_fill_gradient2(expression('weights'),
                           low = "#8b0000",
                           mid = "white",
                           high = "navy",
                           midpoint = 0,
                           breaks = c(-0.001, -0.0005, 0.0005, 0.001, 0.002),
                           labels = c("-0.001", "-0.0005", "0.0005", "0.001", "0.002"),
                           na.value = "grey",
                           limits = c(-0.002, 0.004)) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            line = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text.align = 0.75,
            legend.key.width = unit(100, "points"),
            panel.grid.major = element_line(colour = "transparent"))
    # ggplot() +
    #   xlim(-125, -65) +
    #   ylim(25, 50) +
    #   geom_sf(data = mapdat[mapdat$pm25_binary == 0, ], aes(fill = wts_pooled), color=NA, size = 0.005) +
    #   scale_fill_gradient2(expression('weights'),
    #                        low = "white",
    #                        mid = "white",
    #                        high = "#8b0000",
    #                        midpoint = 0,
    #                        breaks = c(-0.001, -0.0005, 0.002, 0.003, 0.004),
    #                        labels = c("-0.001", "-0.0005", "0.002", "0.003", "0.004"),
    #                        na.value = "grey",
    #                        limits = c(-0.002, 0.004)) +
    #   new_scale_fill() +
    #   geom_sf(data = mapdat[mapdat$pm25_binary == 1, ], aes(fill = wts_pooled), color=NA, size = 0.005) +
    #   scale_fill_gradient2(expression('weights'),
    #                        low = "white",
    #                        mid = "white",
    #                        high = "navy",
    #                        midpoint = 0,
    #                        breaks = c(-0.001, -0.0005, 0.002, 0.003, 0.004),
    #                        labels = c("-0.001", "-0.0005", "0.002", "0.003", "0.004"),
    #                        na.value = "grey",
    #                        limits = c(-0.002, 0.004)) +
    #   geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
    #   theme_minimal() +
    #   theme(
    #     axis.text.x = element_blank(),
    #     axis.text.y = element_blank(),
    #     axis.ticks = element_blank(),
    #     line = element_blank(),
    #     axis.title = element_blank(),
    #     legend.position = "bottom",
    #     legend.direction = "horizontal",
    #     legend.text.align = 0.75,
    #     legend.key.width = unit(100, "points"),
    #     panel.grid.major = element_line(colour = "transparent")
    #   ) #+
      #guides(fill = guide_colourbar(barwidth = 1))
  })
  
  output$binaryPlot <- renderPlot({
      #mapdat$wts_re_treated = ifelse(mapdat$pm25_binary ==1, wts_re(), NA)
    ggplot(mapdat) +
      xlim(-125, -65) +
      ylim(25, 50) +
      geom_sf(aes(fill = factor(pm25_binary)), color = NA, size = 0.005) +  # Use factor(Z2) to treat Z2 as a categorical variable
      scale_fill_manual(
        name = "Treatment",
        values = c("0" = "lightblue", "1" = "red"),  # Specify colors for each level
        labels = c("0" = "PM < 9", "1" = "PM > 9")  # Optional: Provide labels for the legend
      ) +
      geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text.align = 0.75,
        legend.key.width = unit(100, "points"),
        panel.grid.major = element_line(colour = "transparent"),
        legend.title = element_text(size = 20),  # Increase the size of the legend title
        legend.text = element_text(size = 12),   # Increase the size of the legend text
        legend.key.size = unit(2, "lines")       # Increase the size of the legend keys
      )
    
    
    # mapdat$wts_re_treated = ifelse(mapdat$pm25_binary ==1, 1, 0)
    # ggplot(mapdat) +
    #   xlim(-125, -65) +
    #   ylim(25, 50) +
    #   geom_sf(aes(fill = wts_re_treated), color=NA, size = 0.005) +
    #   geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
    #   theme_minimal() +
    #   scale_fill_gradient2(expression('weights'),
    #                        low = "#8b0000",
    #                        mid = "white",
    #                        high = "navy",
    #                        midpoint = 0.5,
    #                        breaks = c(0, 0.25, 0.5, 0.75, 1),
    #                        labels = c("untreated", "", "", "", "treated"),
    #                        na.value = "grey",
    #                        limits = c(-2, 3)
    #                        ) +
    #   theme(#plot.title = element_text(size = 24 * 2,hjust = 0.5),
    #         axis.text.x = element_blank(),
    #         axis.text.y = element_blank(),
    #         axis.ticks = element_blank(),
    #         line = element_blank(),
    #         axis.title = element_blank(),
    #         legend.position = "bottom",
    #         legend.direction = "horizontal",
    #         #legend.text = element_text(angle = 60,  size = 20 * 2),
    #         legend.text.align = 0.75,
    #         #legend.title = element_text(size = 24 * 2),
    #         legend.key.width = unit(100, "points"),
    #         panel.grid.major = element_line(colour = "transparent"))
  })
  
  output$multiPlot <- renderPlot({
    ggplot(mapdat) +
      xlim(-125, -65) +
      ylim(25, 50) +
      geom_sf(aes(fill = factor(pm25_multi)), color = NA, size = 0.005) +  # Use factor(Z2) to treat Z2 as a categorical variable
      scale_fill_manual(
        name = "Treatment",
        values = c("1" = "lightblue", "2" = "green", "3" = "yellow", "4" = "orange", "5" = "red"),  # Specify colors for each level
        labels = c("1" = "PM < 6", "2" = "6 < PM < 9", "3" = "9 < PM < 10", "4" = "10 < PM < 11", "5" = "PM > 11")  # Provide labels for the legend
      ) +
      geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text.align = 0.75,
        legend.key.width = unit(50, "points"),
        panel.grid.major = element_line(colour = "transparent"),
        legend.title = element_text(size = 20),  # Increase the size of the legend title
        legend.text = element_text(size = 12),   # Increase the size of the legend text
        legend.key.size = unit(2, "lines")       # Increase the size of the legend keys
      )
  })
  
  output$carPlot <- renderPlot({
    ggplot(mapdat) +
      xlim(-125, -65) +
      ylim(25, 50) +
      geom_sf(aes(fill = carweights), color=NA, size = 0.005) +
      geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
      theme_minimal() +
      scale_fill_gradient2(expression('weights'),
                           low = "#8b0000",
                           mid = "white",
                           high = "navy",
                           midpoint = 0,
                           breaks = c(-0.001, -0.0005, 0.0005, 0.001, 0.002),
                           labels = c("-0.001", "-0.0005", "0.0005", "0.001", "0.002"),
                           na.value = "grey",
                           limits = c(-0.002, 0.005)) +
      theme(#plot.title = element_text(size = 24 * 2,hjust = 0.5),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            line = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal",
            #legend.text = element_text(angle = 60,  size = 20 * 2),
            legend.text.align = 0.75,
            #legend.title = element_text(size = 24 * 2),
            legend.key.width = unit(100, "points"),
            panel.grid.major = element_line(colour = "transparent"))
    # ggplot() +
    #   xlim(-125, -65) +
    #   ylim(25, 50) +
    #   geom_sf(data = mapdat[mapdat$pm25_binary == 0, ], aes(fill = carweights), color=NA, size = 0.005) +
    #   scale_fill_gradient2(expression('weights'),
    #                        low = "white",
    #                        mid = "white",
    #                        high = "#8b0000",
    #                        midpoint = 0,
    #                        breaks = c(-0.001, -0.0005, 0.002, 0.003, 0.004),
    #                        labels = c("-0.001", "-0.0005", "0.002", "0.003", "0.004"),
    #                        na.value = "grey",
    #                        limits = c(-0.0035, 0.01)) +
    #   new_scale_fill() +
    #   geom_sf(data = mapdat[mapdat$pm25_binary == 1, ], aes(fill = carweights), color=NA, size = 0.005) +
    #   scale_fill_gradient2(expression('weights'),
    #                        low = "white",
    #                        mid = "white",
    #                        high = "navy",
    #                        midpoint = 0,
    #                        breaks = c(-0.001, -0.0005, 0.002, 0.003, 0.004),
    #                        labels = c("-0.001", "-0.0005", "0.002", "0.003", "0.004"),
    #                        na.value = "grey",
    #                        limits = c(-0.0035, 0.01)) +
    #   geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
    #   theme_minimal() +
    #   theme(
    #     axis.text.x = element_blank(),
    #     axis.text.y = element_blank(),
    #     axis.ticks = element_blank(),
    #     line = element_blank(),
    #     axis.title = element_blank(),
    #     legend.position = "bottom",
    #     legend.direction = "horizontal",
    #     legend.text.align = 0.75,
    #     legend.key.width = unit(100, "points"),
    #     panel.grid.major = element_line(colour = "transparent")
    #   )
  })
  
  output$carPlot_binary <- renderPlot({
    ggplot(mapdat) +
      xlim(-125, -65) +
      ylim(25, 50) +
      geom_sf(aes(fill = carweights_binary), color=NA, size = 0.005) +
      geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
      theme_minimal() +
      scale_fill_gradient2(expression('weights'),
                           low = "#8b0000",
                           mid = "white",
                           high = "navy",
                           midpoint = 0,
                           breaks = c(-0.001, -0.0005, 0.0005, 0.001, 0.002),
                           labels = c("-0.001", "-0.0005", "0.0005", "0.001", "0.002"),
                           na.value = "grey",
                           limits = c(min(mapdat$carweights_binary), max(mapdat$carweights_binary))) +
      theme(#plot.title = element_text(size = 24 * 2,hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        #legend.text = element_text(angle = 60,  size = 20 * 2),
        legend.text.align = 0.75,
        #legend.title = element_text(size = 24 * 2),
        legend.key.width = unit(100, "points"),
        panel.grid.major = element_line(colour = "transparent"))
    
  })
   
  output$sarPlot <- renderPlot({
    ggplot(mapdat) +
      xlim(-125, -65) +
      ylim(25, 50) +
      geom_sf(aes(fill = sarweights), color=NA, size = 0.005) +
      geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
      theme_minimal() +
      scale_fill_gradient2(expression('weights'),
                           low = "#8b0000",
                           mid = "white",
                           high = "navy",
                           midpoint = 0,
                           breaks = c(-0.002, -0.001, 0.001, 0.002, 0.004),
                           labels = c("-0.002", "-0.001", "0.001", "0.002", "0.004"),
                           na.value = "grey",
                           limits = c(-0.004, 0.01)) +
      theme(#plot.title = element_text(size = 24 * 2,hjust = 0.5),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            line = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal",
            #legend.text = element_text(angle = 60,  size = 20 * 2),
            legend.text.align = 0.75,
            #legend.title = element_text(size = 24 * 2),
            legend.key.width = unit(100, "points"),
            panel.grid.major = element_line(colour = "transparent"))
    # ggplot() +
    #   xlim(-125, -65) +
    #   ylim(25, 50) +
    #   geom_sf(data = mapdat[mapdat$pm25_binary == 0, ], aes(fill = sarweights), color=NA, size = 0.005) +
    #   scale_fill_gradient2(expression('weights'),
    #                        low = "white",
    #                        mid = "white",
    #                        high = "#8b0000",
    #                        midpoint = 0,
    #                        breaks = c(-0.001, -0.0005, 0.002, 0.003, 0.004),
    #                        labels = c("-0.001", "-0.0005", "0.002", "0.003", "0.004"),
    #                        na.value = "grey",
    #                        limits = c(-0.0035, 0.01)) +
    #   new_scale_fill() +
    #   geom_sf(data = mapdat[mapdat$pm25_binary == 1, ], aes(fill = sarweights), color=NA, size = 0.005) +
    #   scale_fill_gradient2(expression('weights'),
    #                        low = "white",
    #                        mid = "white",
    #                        high = "navy",
    #                        midpoint = 0,
    #                        breaks = c(-0.001, -0.0005, 0.002, 0.003, 0.004),
    #                        labels = c("-0.001", "-0.0005", "0.002", "0.003", "0.004"),
    #                        na.value = "grey",
    #                        limits = c(-0.0035, 0.01)) +
    #   geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
    #   theme_minimal() +
    #   theme(
    #     axis.text.x = element_blank(),
    #     axis.text.y = element_blank(),
    #     axis.ticks = element_blank(),
    #     line = element_blank(),
    #     axis.title = element_blank(),
    #     legend.position = "bottom",
    #     legend.direction = "horizontal",
    #     legend.text.align = 0.75,
    #     legend.key.width = unit(100, "points"),
    #     panel.grid.major = element_line(colour = "transparent")
    #   )
  })
  
  output$sarPlot_binary <- renderPlot({
    ggplot(mapdat) +
      xlim(-125, -65) +
      ylim(25, 50) +
      geom_sf(aes(fill = sarweights_binary), color=NA, size = 0.005) +
      geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
      theme_minimal() +
      scale_fill_gradient2(expression('weights'),
                           low = "#8b0000",
                           mid = "white",
                           high = "navy",
                           midpoint = 0,
                           breaks = c(-0.002, -0.001, 0.001, 0.002, 0.004),
                           labels = c("-0.002", "-0.001", "0.001", "0.002", "0.004"),
                           na.value = "grey",
                           limits = c(-0.004, 0.01)) +
      theme(#plot.title = element_text(size = 24 * 2,hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        #legend.text = element_text(angle = 60,  size = 20 * 2),
        legend.text.align = 0.75,
        #legend.title = element_text(size = 24 * 2),
        legend.key.width = unit(100, "points"),
        panel.grid.major = element_line(colour = "transparent"))
  })
  
  output$balance <- renderTable({
    balance = data.frame('Covariate' = colnames(Xstate))
    balance$unadjusted = abs(colMeans(Xstate[mapdat$pm25_binary == 1,])-colMeans(Xstate[mapdat$pm25_binary == 0,]))

    # Calculate the weighted mean of covariates in each treatment group using the weights wts_re
    Xtreated = Xstate[mapdat$pm25_binary == 1,]
    wtspooled_treated = mapdat$wts_pooled[mapdat$pm25_binary == 1]
    wtsre_treated = wts_re()[mapdat$pm25_binary == 1]
    wtsfe_treated = mapdat$wts_fe[mapdat$pm25_binary == 1]
    wtscar_treated = mapdat$carweights[mapdat$pm25_binary == 1]
    wtsar_treated = mapdat$sarweights[mapdat$pm25_binary == 1]
    
    Xcontrol = Xstate[mapdat$pm25_binary == 0,]
    wtspooled_control = mapdat$wts_pooled[mapdat$pm25_binary == 0]
    wtsre_control = wts_re()[mapdat$pm25_binary == 0]
    wtsfe_control = mapdat$wts_fe[mapdat$pm25_binary == 0]
    wtscar_control = mapdat$carweights[mapdat$pm25_binary == 0]
    wtsar_control = mapdat$sarweights[mapdat$pm25_binary == 0]
    
    balance$pooled = abs(t(Xtreated) %*% wtspooled_treated - t(Xcontrol) %*% wtspooled_control)
    balance$randeff = abs(t(Xtreated) %*% wtsre_treated - t(Xcontrol) %*% wtsre_control)
    balance$fixedeff = abs(t(Xtreated) %*% wtsfe_treated - t(Xcontrol) %*% wtsfe_control)
    balance$car = abs(t(Xtreated) %*% wtscar_treated - t(Xcontrol) %*% wtscar_control)
    balance$sar = abs(t(Xtreated) %*% wtsar_treated - t(Xcontrol) %*% wtsar_control)
    
    balance = balance[-1,]
    vartreated = apply(Xtreated, 2, var)
    varcontrol = apply(Xcontrol, 2, var)
    balance$pooled = round(balance$pooled/sqrt((vartreated[-1] + varcontrol[-1])/2),5)
    balance$randeff = round(balance$randeff/sqrt((vartreated[-1] + varcontrol[-1])/2),5)
    balance$fixedeff = round(balance$fixedeff/sqrt((vartreated[-1] + varcontrol[-1])/2),5)
    balance$unadjusted = round(balance$unadjusted/sqrt((vartreated[-1] + varcontrol[-1])/2),5)
    balance$car = round(balance$car/sqrt((vartreated[-1] + varcontrol[-1])/2),5)
    balance$sar = round(balance$sar/sqrt((vartreated[-1] + varcontrol[-1])/2),5)
    
    global$balance = balance
    balance[1:(ncol(X)-1),]
  })

  output$statebalance <- renderTable({
    global$balance[ncol(X):nrow(global$balance),]
  })
  
  output$estimates <- renderTable({
    # initialize empty data frame of 8 columns
    estimates = data.frame('Unadjusted' = rep(NA,2),
                           'Pooled' = rep(NA,2),
                           'RE' = rep(NA,2),
                           'FE' = rep(NA,2),
                           'CAR' = rep(NA,2),
                           'SAR' = rep(NA,2),
                           'CARBinary' = rep(NA,2),
                           'SARBinary' = rep(NA,2))
    
    estimates$Unadjusted[1] = mean(Y[mapdat$pm25_binary == 1])-mean(Y[mapdat$pm25_binary == 0])
    estimates$Pooled[1] = sum(Y[mapdat$pm25_binary == 1]*mapdat$wts_pooled[mapdat$pm25_binary == 1]) - sum(Y[mapdat$pm25_binary == 0]*mapdat$wts_pooled[mapdat$pm25_binary == 0])
    estimates$RE[1] = sum(Y[mapdat$pm25_binary == 1]*wts_re()[mapdat$pm25_binary == 1]) - sum(Y[mapdat$pm25_binary == 0]*wts_re()[mapdat$pm25_binary == 0])
    estimates$FE[1] = sum(Y[mapdat$pm25_binary == 1]*mapdat$wts_fe[mapdat$pm25_binary == 1]) - sum(Y[mapdat$pm25_binary == 0]*mapdat$wts_fe[mapdat$pm25_binary == 0])
    estimates$CAR[1] = sum(Y[mapdat$pm25_binary == 1]*mapdat$carweights[mapdat$pm25_binary == 1]) - sum(Y[mapdat$pm25_binary == 0]*mapdat$carweights[mapdat$pm25_binary == 0])
    estimates$SAR[1] = sum(Y[mapdat$pm25_binary == 1]*mapdat$sarweights[mapdat$pm25_binary == 1]) - sum(Y[mapdat$pm25_binary == 0]*mapdat$sarweights[mapdat$pm25_binary == 0])
    estimates$CARBinary[1] = sum(Y[mapdat$pm25_binary == 1]*mapdat$carweights_binary[mapdat$pm25_binary == 1]) - sum(Y[mapdat$pm25_binary == 0]*mapdat$carweights_binary[mapdat$pm25_binary == 0])
    estimates$SARBinary[1] = sum(Y[mapdat$pm25_binary == 1]*mapdat$sarweights_binary[mapdat$pm25_binary == 1]) - sum(Y[mapdat$pm25_binary == 0]*mapdat$sarweights_binary[mapdat$pm25_binary == 0])
    
    mod_pooled = lm(all_all_65up_imputed ~ -1 + pm25_binary + X, data = mapdat)
    estimates$Pooled[2] = mod_pooled$coefficients['pm25_binary']
    mod_fe = lm(all_all_65up_imputed ~ -1 + pm25_binary + X + factor(statefactor), data = mapdat)
    estimates$FE[2] = mod_fe$coefficients['pm25_binary']
    mod_re = glmmTMB(all_all_65up_imputed ~ -1 + pm25_binary + X + (1|statefactor), data = mapdat, 
                     start = list(betad = log(input$sig2eps)/2, theta = log(input$sig2gam)/2),
                     map = list(betad = factor(NA), theta = factor(NA)))
    estimates$RE[2] = summary(mod_re)$coefficients$cond['pm25_binary','Estimate']
    rownames(estimates) = c('Using Weights', 'Using lm, glmmTMB')
    estimates
  }, rownames = TRUE)
  
  output$multiestimates <- renderTable({
    estimates = data.frame('Single Model, all Interactions' = rep(NA,4),
                           'Separate models' = rep(NA,4),
                           'Common Weighting' = rep(NA,4))
    
    Xsub = st_drop_geometry(mapdat[,20:32])
    Xsub = as.matrix(Xsub)
    
    mapdat$Z1 = 1 * (mapdat$pm25_multi == 1)
    mapdat$Z2 = 1 * (mapdat$pm25_multi == 2)
    mapdat$Z3 = 1 * (mapdat$pm25_multi == 3)
    mapdat$Z4 = 1 * (mapdat$pm25_multi == 4)
    mapdat$Z5 = 1 * (mapdat$pm25_multi == 5)
    
    Xsub_centered = scale(Xsub, center = TRUE, scale = FALSE)
    # Solution 1: Fit a single model with interaction terms
    mod_interacted = lm(all_all_65up_imputed ~ (Z2 + Z3 + Z4 + Z5)*Xsub_centered, data = mapdat)
    # paste coefficient and stnadard error for each term as Est(SE)
    est = round(summary(mod_interacted)$coefficients['Z2', c('Estimate', 'Std. Error')],3)
    estimates[1,1] = paste0(est[1], '(', est[2], ')')
    est = round(summary(mod_interacted)$coefficients['Z3', c('Estimate', 'Std. Error')],3)
    estimates[2,1] = paste0(est[1], '(', est[2], ')')
    est = round(summary(mod_interacted)$coefficients['Z4', c('Estimate', 'Std. Error')],3)
    estimates[3,1] = paste0(est[1], '(', est[2], ')')
    est = round(summary(mod_interacted)$coefficients['Z5', c('Estimate', 'Std. Error')],3)
    estimates[4,1] = paste0(est[1], '(', est[2], ')')
    
    # Solution 2: Fit separate models for each treatment arm
    mod_Z2 = lm(all_all_65up_imputed ~ Z2 + Xsub, data = mapdat)
    est = round(summary(mod_Z2)$coefficients['Z2', c('Estimate', 'Std. Error')],3)
    estimates[1,2] = paste0(est[1], '(', est[2], ')')
    mod_Z3 = lm(all_all_65up_imputed ~ Z3 + Xsub, data = mapdat)
    est = round(summary(mod_Z3)$coefficients['Z3', c('Estimate', 'Std. Error')],3)
    estimates[2,2] = paste0(est[1], '(', est[2], ')')
    mod_Z4 = lm(all_all_65up_imputed ~ Z4 + Xsub, data = mapdat)
    est = round(summary(mod_Z4)$coefficients['Z4', c('Estimate', 'Std. Error')],3)
    estimates[3,2] = paste0(est[1], '(', est[2], ')')
    mod_Z5 = lm(all_all_65up_imputed ~ Z5 + Xsub, data = mapdat)
    est = round(summary(mod_Z5)$coefficients['Z5', c('Estimate', 'Std. Error')],3)
    estimates[4,2] = paste0(est[1], '(', est[2], ')')
    
    # Solution 3: Common weighting
    pi1 = mean(mapdat$pm25_multi == 1)
    pi2 = mean(mapdat$pm25_multi == 2)
    pi3 = mean(mapdat$pm25_multi == 3)
    pi4 = mean(mapdat$pm25_multi == 4)
    pi5 = mean(mapdat$pm25_multi == 5)
    pi = c(pi1, pi2, pi3, pi4, pi5)
    # Fit multinomial logistic regression
    propensity = multinom(pm25_multi ~ Xsub, data = mapdat)
    # predict probabilities for each unit using propensity # phatk(Wi)
    prob = predict(propensity, newdata = mapdat, type = 'probs') 
    # Replace any zero entries in prob with 0.000001
    prob[prob == 0] = 0.000001
    # compute lamhatinv
    lamhatinv = rep(NA, nrow(mapdat))
    for (i in 1:nrow(mapdat)) {
      lamhatinv[i] = sum(pi*(1-pi)/prob[i,])
    }
    lamhat = 1/lamhatinv
    estimates[1,3] = round((1/sum(lamhat*mapdat$Z2/prob[,2]))*sum(lamhat*mapdat$Z2*mapdat$all_all_65up_imputed/prob[,2]) - 
      (1/sum(lamhat*mapdat$Z1/prob[,1]))*sum(lamhat*mapdat$Z1*mapdat$all_all_65up_imputed/prob[,1]),3)
    estimates[2,3] = round((1/sum(lamhat*mapdat$Z3/prob[,3]))*sum(lamhat*mapdat$Z3*mapdat$all_all_65up_imputed/prob[,3]) - 
      (1/sum(lamhat*mapdat$Z1/prob[,1]))*sum(lamhat*mapdat$Z1*mapdat$all_all_65up_imputed/prob[,1]),3)
    estimates[3,3] = round((1/sum(lamhat*mapdat$Z4/prob[,4]))*sum(lamhat*mapdat$Z4*mapdat$all_all_65up_imputed/prob[,4]) - 
      (1/sum(lamhat*mapdat$Z1/prob[,1]))*sum(lamhat*mapdat$Z1*mapdat$all_all_65up_imputed/prob[,1]),3)
    estimates[4,3] = round((1/sum(lamhat*mapdat$Z5/prob[,5]))*sum(lamhat*mapdat$Z5*mapdat$all_all_65up_imputed/prob[,5]) - 
      (1/sum(lamhat*mapdat$Z1/prob[,1]))*sum(lamhat*mapdat$Z1*mapdat$all_all_65up_imputed/prob[,1]),3)
    
    
    rownames(estimates) = c('Treatment 2', 'Treatment 3', 'Treatment 4', 'Treatment 5')
    estimates
  }, rownames = TRUE)
  
  output$dispersion <- renderTable({
    # Initialize empty data frame with column names pooled, randeff, fixedeff, car, sar
    dispersion = data.frame('Method' = c('Pooled', 'Random Effects', 'Fixed Effects', 'CAR', 'SAR'))
    dispersion$dispersion = c(var(mapdat$wts_pooled), var(wts_re()), var(mapdat$wts_fe), 
                              var(mapdat$carweights), var(mapdat$sarweights))
    # Convert dispersion column to scientific notation
    dispersion$dispersion = format(dispersion$dispersion, scientific = TRUE)
    # Now compute 
    dispersion$ESS = c(sum(abs(mapdat$wts_pooled))^2/sum(mapdat$wts_pooled^2),
                       sum(abs(wts_re()))^2/sum(wts_re()^2),
                       sum(abs(mapdat$wts_fe))^2/sum(mapdat$wts_fe^2),
                       sum(abs(mapdat$carweights))^2/sum(mapdat$carweights^2),
                       sum(abs(mapdat$sarweights))^2/sum(mapdat$sarweights^2))
    dispersion
  })
  
  output$weighthists <- renderPlot({
    data <- data.frame(
      group = rep(c("Pooled", "Random Effects", "Fixed Effects", "CAR", "SAR"), each = nrow(mapdat)),
      value = c(mapdat$wts_pooled, wts_re(), mapdat$wts_fe, mapdat$carweights, mapdat$sarweights)
    )
    data$group <- factor(data$group, levels = c("Pooled", "Random Effects", "Fixed Effects", "CAR", "SAR"))
    
    x_limits <- range(data$value, na.rm = TRUE)
    
    # Create multiple histograms using facets
    ggplot(data, aes(x = value)) +
      geom_histogram(fill = "lightblue") +
      labs(x = "weights", y = "Frequency") +
      facet_wrap(~ group, nrow = 1) +
      theme_minimal() + 
      coord_cartesian(xlim = x_limits)
  })
  # 
  # output$empiricalBayes <- renderTable({
  #   Y = mapdat$cms_mortality_pct
  #   eb = data.frame('State' = colnames(state.dummy))
  #   eb$empiricalbayes = empirical_bayes(Y = Y, 
  #                                       X=X, 
  #                                       Z=Z, 
  #                                       clusters = clusters, 
  #                                       sig2gam = input$sig2gam, 
  #                                       sig2eps = input$sig2eps)
  #   eb$empiricalbayes =  eb$empiricalbayes*10000
  #   eb
  # })
  
  # Create a table with columns corresponding to 
  # that displays the within-state balance ASMDs for each covariate
  output$balancebystate <- renderTable({
    balancebystate = data.frame('Covariate' = colnames(X))
    # for each state, calculate the ASMD for each covariate
    for (i in 1:length(colnames(state.dummy))) {
      st = colnames(state.dummy)[i]
      Xtreated = matrix(X[mapdat$pm25_binary == 1 & mapdat$STATE_CODE == st,], ncol = ncol(X))
      
      wtsre_treated = wts_re()[mapdat$pm25_binary == 1 & mapdat$STATE_CODE == st]
      if (nrow(Xtreated) == 0){
        treated_contribution = rep(0,ncol(X))
        vartreated = rep(0,ncol(X))
      }
      else{
        treated_contribution = t(Xtreated) %*% wtsre_treated
        if (nrow(Xtreated) == 1){
          vartreated = rep(0,ncol(X))
        }
        else{
          vartreated = apply(Xtreated, 2, var)
        }
      }
      Xcontrol = matrix(X[mapdat$pm25_binary == 0 & mapdat$STATE_CODE == st,], ncol = ncol(X))
      wtsre_control = wts_re()[mapdat$pm25_binary == 0 & mapdat$STATE_CODE == st]
      if (nrow(Xcontrol) == 0){
        control_contribution = rep(0,ncol(X))
        varcontrol = rep(0,ncol(X))
      }
      else{
        control_contribution = t(Xcontrol) %*% wtsre_control
        if (nrow(Xcontrol) == 1){
          varcontrol = rep(0,ncol(X))
        }
        else{
          varcontrol = apply(Xcontrol, 2, var)
        }
      }
      
      balancebystate[[st]] = abs(treated_contribution - control_contribution)
      
      # 
      # )
      balancebystate[[st]] = round(balancebystate[[st]]/sqrt((vartreated + varcontrol)/2),5)
    }
    # Remove DC
    balancebystate = balancebystate[, -which(colnames(balancebystate) == 'DC')]
    # Remove Intercept
    balancebystate = balancebystate[-1,]
    balancebystate
  })
  
  output$outcomeplot <- renderPlot({
    ggplot(mapdat) +
      xlim(-125, -65) +
      ylim(25, 50) +
      geom_sf(aes(fill = all_all_65up), color=NA, size = 0.005) +
      geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
      theme_minimal() +
      scale_fill_gradient2(expression('weights'),
                           low = "#8b0000",
                           mid = "white",
                           high = "navy",
                           midpoint = median(mapdat$all_all_65up, na.rm = TRUE),
                           breaks = quantile(mapdat$all_all_65up, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = T),
                           labels = as.character(round(quantile(mapdat$all_all_65up, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = T), 3)),
                           na.value = "grey",
                           limits = c(min(mapdat$all_all_65up), max(mapdat$all_all_65up))) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            line = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text.align = 0.75,
            legend.key.width = unit(100, "points"),
            panel.grid.major = element_line(colour = "transparent"))
  })
  
  output$outcome_imputed_plot <- renderPlot({
    ggplot(mapdat) +
      xlim(-125, -65) +
      ylim(25, 50) +
      geom_sf(aes(fill = all_all_65up_imputed), color=NA, size = 0.005) +
      geom_sf(data = stateboundaries, fill = NA, color = "black", size = 3) +
      theme_minimal() +
      scale_fill_gradient2(expression('weights'),
                           low = "#8b0000",
                           mid = "white",
                           high = "navy",
                           midpoint = median(mapdat$all_all_65up_imputed, na.rm = TRUE),
                           breaks = quantile(mapdat$all_all_65up_imputed, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = T),
                           labels = as.character(round(quantile(mapdat$all_all_65up_imputed, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = T), 3)),
                           na.value = "grey",
                           limits = c(min(mapdat$all_all_65up_imputed), max(mapdat$all_all_65up_imputed))) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            line = element_blank(),
            axis.title = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.text.align = 0.75,
            legend.key.width = unit(100, "points"),
            panel.grid.major = element_line(colour = "transparent"))
  })
  # output$balancebystate_QX <- renderTable({
  #   # siginvsq = list()
  #   # for (k in 1:length(unique(clusters))){
  #   #   sigk = input$sig2eps*diag(1, nks[k]) + input$sig2gam*matrix(1, nks[k], nks[k])
  #   #   siginvsq[[k]] = sqrtm(solve(sigk))
  #   # }
  #   # Q = bdiag(siginvsq)
  #   load('sqre.RData')
  #   Q = sqre
  #   QX = Q %*% X
  #   QZ = Q %*% Z
  #   balancebystate_QX = data.frame('Covariate' = colnames(QX))
  #   # for each state, calculate the ASMD for each covariate
  #   for (i in 1:length(colnames(state.dummy))) {
  #     st = colnames(state.dummy)[i]
  #     QZstate = QZ[mapdat$STATE_CODE == st]
  #     Zstate = Z[mapdat$STATE_CODE == st]
  #     QXstate = matrix(QX[mapdat$STATE_CODE == st,], ncol = ncol(QX))
  #     QXtreated = matrix(QXstate[Zstate == 1,], ncol = ncol(QX))
  #     wtsre_treated = wts_re()[mapdat$pm25_binary == 1 & mapdat$STATE_CODE == st]
  #     if (nrow(QXtreated) == 0){
  #       treated_contribution = rep(0,ncol(QX))
  #       vartreated = rep(0,ncol(QX))
  #     }
  #     else{
  #       treated_contribution = t(QXtreated) %*% wtsre_treated
  #       if (nrow(QXtreated) == 1){
  #         vartreated = rep(0,ncol(QX))
  #       }
  #       else{
  #         vartreated = apply(QXtreated, 2, var)
  #       }
  #     }
  #     QXcontrol = matrix(QXstate[Zstate == 0,], ncol = ncol(QX))
  #     wtsre_control = wts_re()[mapdat$pm25_binary == 0 & mapdat$STATE_CODE == st]
  #     if (nrow(QXcontrol) == 0){
  #       control_contribution = rep(0,ncol(QX))
  #       varcontrol = rep(0,ncol(QX))
  #     }
  #     else{
  #       control_contribution = t(QXcontrol) %*% wtsre_control
  #       if (nrow(QXcontrol) == 1){
  #         varcontrol = rep(0,ncol(QX))
  #       }
  #       else{
  #         varcontrol = apply(QXcontrol, 2, var)
  #       }
  #     }
  #     QZpos = max(QZstate[Zstate == 1], 0, na.rm = T)
  #     QZneg = min(QZstate[Zstate == 0], 0, na.rm = T)
  #     balancebystate_QX[[st]] = abs(QZpos[1]*treated_contribution - abs(QZneg[1])*control_contribution)
  #     
  #     # 
  #     # )
  #     balancebystate_QX[[st]] = round(balancebystate_QX[[st]]/sqrt((vartreated + varcontrol)/2),5)
  #   }
  #   # Remove DC
  #   balancebystate_QX = balancebystate_QX[, -which(colnames(balancebystate_QX) == 'DC')]
  #   # Remove Intercept
  #   balancebystate_QX = balancebystate_QX[-1,]
  #   balancebystate_QX
  # })
}

shinyApp(ui = ui, server = server)