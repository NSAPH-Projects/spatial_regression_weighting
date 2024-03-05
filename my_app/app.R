library(shiny)

source('/Users/sophie/Documents/implied_weights/impliedweights_randeffs.R')

## Import data
study = read.csv('/Users/sophie/Documents/SpatialConf/archived/Study_dataset_2010.csv')

study$FIPS = str_pad(study$FIPS, 5, pad = '0')
study$pm25_binary = ifelse(study$qd_mean_pm25 > 9, 1, 0)

# Create design matrix. 
study$statefactor = as.numeric(as.factor(study$STATE_CODE))
#studyexclude = which(study$STATE_CODE %in% c("CA", "OR", "WA", "ID", "NV", "UT", "AZ", "CO", "NM", "WY", "MT"))
#study = study[-studyexclude,]

study = mutate(study,
               STATEFP10 = str_sub(FIPS, 1, 2),
               COUNTYFP10 = str_sub(FIPS, 3, 5))

# Import adjacency matrix created from spacebench script
#adj = read.csv("/Users/sophie/Documents/SpatialConf/archived/adjacency_matrix.csv", header = F)

# column bind adj to study
#study = cbind(study, adj[-studyexclude,-studyexclude])

states = st_read("/Users/sophie/Documents/SpatialConf/archived/tl_2010_us_county10/tl_2010_us_county10.shp")
stateboundaries = st_read("/Users/sophie/Downloads/tl_2010_us_state10/tl_2010_us_state10.shp")

# Omit hispanic (collinear with other races)
mapdat = right_join(states, study, by = c("STATEFP10", "COUNTYFP10"))
# Order mapdat by statefactor
mapdat = mapdat[order(mapdat$statefactor),]
clusters = mapdat$statefactor

# Extract A as the last 2695 columns of mapdat
#A = st_drop_geometry(mapdat[,(ncol(mapdat)-2695):(ncol(mapdat)-1)])
#A = as.matrix(A)
X = st_drop_geometry(mapdat[,c(23:29,33:34)])
X = cbind(1, X)
X = as.matrix(X)
Z = mapdat$pm25_binary
state.dummy = model.matrix(~STATE_CODE - 1, data = mapdat)
colnames(state.dummy)=gsub("STATE_CODE","",as.character(colnames(state.dummy)))
Xstate = cbind(1, st_drop_geometry(mapdat[,c(23:29,33:34)]), state.dummy[,-1])
Xstate = as.matrix(Xstate)

mapdat$wts_fe = fixedeffects(
  X = Xstate,
  Z = Z
)
load("../carweights07.rdata")
load("../sarweights07.rdata")
mapdat$carweights = carweights
mapdat$sarweights = sarweights

p = ncol(X)
K = length(unique(clusters))
Xks = lapply(1:K, function(k) X[clusters == k,])
Zks = lapply(1:K, function(k) Z[clusters == k])
nks = sapply(1:K, function(k) sum(clusters == k))
mks = sapply(1:K, function(k) sum(Z[clusters == k]))

Xbarks = lapply(1:K, function(k) {
  Xk = matrix(X[clusters == k,], ncol = p, nrow = nks[k])
  return(matrix(colMeans(Xk), ncol = p, nrow = 1))
}
)

Xtbarks = lapply(1:K, function(k) {
  if (mks[k] > 0){
    if (nks[k] == 1){
      return(Xks[[k]])
    }
    else{
      Xktreated = matrix(Xks[[k]][Zks[[k]] == 1,], nrow = mks[k], ncol = p)
      return(matrix(colMeans(Xktreated), ncol = p, nrow = 1))
    }
  }
  else{
    return(matrix(0, ncol = p, nrow = 1))
  }
}
)


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel(h1("The Trinity Explained: Implied weights of Fixed Effect, Random Effect and Autoregressive Spatial models", align = "center")),
  
  sidebarPanel(
    # Input: Slider for sig2gam ----
    sliderInput(inputId = "sig2gam",
                label = "sig2gam:",
                min = 0,
                max = 10,
                value = 0.5),
    # Input: Slider for sig2eps ----
    sliderInput(inputId = "sig2eps",
                label = "sig2eps:",
                min = 0.00001,
                max = 10,
                value = 7),
    # Input: slider for phi ----
    # sliderInput(inputId = "phi",
    #             label = "phi:",
    #             min = 0,
    #             max = 0.999,
    #             value = 0.5),
    checkboxInput("fe_weights", "Show Fixed Effects implied weights", TRUE)
  ),
  
  tabsetPanel(
    tabPanel("RE/FE",
             fluidRow(
               column(6, #style='padding:5px;',      
                      h4("Random Effects", align = "center"),
                      plotOutput(outputId = "rePlot")),
               column(6, #style='padding:5px;',
                      h4("Fixed Effects", align = "center"),
                      plotOutput(outputId = "fePlot"))
             )
             ),
    tabPanel("Treated/Control",
             fluidRow(
               column(6,
                      h4("Treated Weights", align = "center"),
                      plotOutput(outputId = "retreatedPlot")),
               column(6,
                      h4("Control Weights", align = "center"),
                      plotOutput(outputId = "recontrolPlot"))
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
    tabPanel("Balance",
             fluidRow(
               column(5,
                      align = "center",
                      h4("Covariate ASMD", align = "center"),
                      tableOutput("balance")),
               column(5,
                      align = "center",
                      h4("State-Indicator ASMD", align = "center"),
                      tableOutput("statebalance")),
               column(2,
                      align = "center",
                      h4("Empirical Bayes (x 10^4)", align = "center"),
                      tableOutput("empiricalBayes"))
             )
            )
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
  
  output$fePlot <- renderPlot({
    if (input$fe_weights) {
      

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
    }
  })
  
  output$retreatedPlot <- renderPlot({
      mapdat$wts_re_treated = ifelse(mapdat$pm25_binary ==1, wts_re(), NA)
      ggplot(mapdat) +
        xlim(-125, -65) + 
        ylim(25, 50) +
        geom_sf(aes(fill = wts_re_treated), color=NA, size = 0.005) +
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
  
  output$recontrolPlot <- renderPlot({
    mapdat$wts_re_control = ifelse(mapdat$pm25_binary==0, wts_re(), NA)
    ggplot(mapdat) +
      xlim(-125, -65) + 
      ylim(25, 50) +
      geom_sf(aes(fill = wts_re_control), color=NA, size = 0.005) +
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
  })
  
  output$balance <- renderTable({
    balance = data.frame('Covariate' = colnames(Xstate))
    balance$unadjusted = abs(colMeans(Xstate[mapdat$pm25_binary == 1,])-colMeans(Xstate[mapdat$pm25_binary == 0,]))

    # Calculate the weighted mean of covariates in each treatment group using the weights wts_re
    Xtreated = Xstate[mapdat$pm25_binary == 1,]
    wtsfe_treated = mapdat$wts_fe[mapdat$pm25_binary == 1]
    wtsre_treated = wts_re()[mapdat$pm25_binary == 1]
    wtscar_treated = mapdat$carweights[mapdat$pm25_binary == 1]
    wtsar_treated = mapdat$sarweights[mapdat$pm25_binary == 1]
    
    Xcontrol = Xstate[mapdat$pm25_binary == 0,]
    wtsfe_control = mapdat$wts_fe[mapdat$pm25_binary == 0]
    wtsre_control = wts_re()[mapdat$pm25_binary == 0]
    wtscar_control = mapdat$carweights[mapdat$pm25_binary == 0]
    wtsar_control = mapdat$sarweights[mapdat$pm25_binary == 0]
    
    balance$randeff = abs(t(Xtreated) %*% wtsre_treated - t(Xcontrol) %*% wtsre_control)
    balance$fixedeff = abs(t(Xtreated) %*% wtsfe_treated - t(Xcontrol) %*% wtsfe_control)
    balance$car = abs(t(Xtreated) %*% wtscar_treated - t(Xcontrol) %*% wtscar_control)
    balance$sar = abs(t(Xtreated) %*% wtsar_treated - t(Xcontrol) %*% wtsar_control)
    balance = balance[-1,]
    vartreated = apply(Xtreated, 2, var)
    varcontrol = apply(Xcontrol, 2, var)
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
  # 
  output$empiricalBayes <- renderTable({
    Y = mapdat$cms_mortality_pct
    eb = data.frame('State' = colnames(state.dummy))
    eb$empiricalbayes = empirical_bayes(Y = Y, 
                                        X=X, 
                                        Z=Z, 
                                        clusters = clusters, 
                                        sig2gam = input$sig2gam, 
                                        sig2eps = input$sig2eps)
    eb$empiricalbayes =  eb$empiricalbayes*10000
    eb
  })
  
}

shinyApp(ui = ui, server = server)