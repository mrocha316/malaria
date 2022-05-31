library(shiny)



####I want to add dynamic tables for alpha0-n and psi1-n
####Look at here for reference
#### https://gist.github.com/christopherlovell/b7ecdf8b0aa82c20fa46
#### https://stackoverflow.com/questions/19130455/create-dynamic-number-of-input-elements-with-r-shiny



inline = function (x) {
  tags$div(style="display:inline-block;", x)
}

# Define UI for miles per gallon app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("n Treatment Bins Malaria Model"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      
      # Input: Slider for variable to choose number of bins ----
      sliderInput("bins",
                  "Number of Treatment Bins:",
                  min=0,
                  max=10,
                  value=0),
      
      # Input: Vector for alpha values
      textInput('vec1', "Treatment rates to different bins (separete by commas)",
                "0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05"),
      
      # Input: Vector for psi values
      textInput('vec2', "Transmission Reduction Factor for treatment groups (sepate by commas)",
                "0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50"),
      
      # Input: Vector for starting conditions
      textInput('vec3', "Starting Proportions Hs, Hi, Hr, Ms, Mi, T_1, ..., T_n (comma sepated)",
                "0.8, 0.1, 0.1, 0.8, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0"),
      
      #Input with side-by-side parameters

      # inline(textInput("popH", "Human Starting Population:", value=35)),
      # inline(textInput("popM", "Mosquito Starting Population:", value=50)),
      
      div(
        style = "display: flex; flex-wrap: wrap;",
        div(
          style = "flex: 1;",
          sliderInput("poph", 
                      "Human Starting Population:", 
                      min=0, max= 100000000, step=5000000,
                      value=35000000)
        ),
        div(
          style = "flex: 1;",
          sliderInput("popm", 
                      "Mosquito Starting Population:", 
                      min=0, max= 200000000,step=10000000,
                      value=105000000)
        ),
        
      ),
      
      # Input: Numbers for parameters
      # sliderInput("poph", 
      #              "Human Starting Population:", 
      #              min=0, max= 100000000, step=5000000,
      #              value=35000000),
      # 
      # sliderInput("popm", 
      #              "Mosquito Starting Population:", 
      #              min=0, max= 200000000,step=10000000,
      #              value=105000000),
      
      div(
        style = "display: flex; flex-wrap: wrap;",
        div(
          style = "flex: 1;",
          numericInput("MuH", 
                       "Human Natural Death Rate:", 
                       min=0, max=1, step=0.0000001,
                       0.0000443)
        ),
        div(
          style = "flex: 1;",
          numericInput("MuM", 
                       "Mosquito Natural Death Rate:", 
                       min=0, max=1, step=0.0000001,
                       1/14)
        ),
        
      ),
      
      div(
        style = "display: flex; flex-wrap: wrap;",
        div(
          style = "flex: 1;",
          numericInput("LamH", 
                       "Human Births:", 
                       1550.5)
        ),
        div(
          style = "flex: 1;",
          numericInput("LamM", 
                       "Mosquito Births:", 
                       7500000)
        ),
        
      ),
      
      div(
        style = "display: flex; flex-wrap: wrap;",
        div(
          style = "flex: 1;",
          numericInput("betaH", 
                        "Rate of moquito to human infection:", 
                        min=0, max=1, step=0.0000001,
                        0.02)
        ),
        div(
          style = "flex: 1;",
          numericInput("betaM", 
                       "Rate of human to mosquito infection:", 
                       min=0, max=1, step=0.0000001,
                       0.01)
        ),
        
      ),
      
            numericInput("omega", 
                   "Rate that humans lose natural immunity:", 
                   min=0, max=1, step=0.0000001,
                   1/28),
      
      numericInput("sigma", 
                   "Rate of natural recovery from infection:", 
                   min=0, max=1, step=0.0000001,
                   1/33),
      
      numericInput("delS", 
                   "Additional disease-induced death rate:", 
                   min=0, max=1, step=0.0000001,
                   0.0000143),
      
      
      numericInput("step", 
                   "Precision of each time step:", 
                   min=0, max=1, step=0.01,
                   0.01),
      
      numericInput("time1",
                   "Total time (days) to model:",
                   min=0, 
                   500)
      
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Result Matrix to test----
      uiOutput('results'),
      
      # Output: Formatted text for caption ----
      h3(textOutput("caption1")),
      
      # Output: Plot of the requested variable against mpg ----
      plotOutput("HumanPlot"),
      
      # Output: Formatted text for caption ----
      h3(textOutput("caption2")),
      
      # Output: Plot of the requested variable against mpg ----
      plotOutput("MosquitoPlot"),
      

    )
  )
)


# Data pre-processing ---- adding the main function
# since this doesn't rely on any user inputs, we can do this once at 
# startup and then use the value throughout the lifetime of the app

fthreeN <- function(first, iterations = 50, n=0, conv = 1*exp(-16), LamH, LamM, betaM, betaH,
                    step, MuH, delS, sigma, omega, MuM, poph, popm, alpha, psi){
  x = matrix(0, nrow=iterations, ncol=5+n)
  x[1,] = first
  if(n==0){
    colnames(x) <- c("Sh", "Ih", "Rh", "Sm", "Im")
    gam <- function(x){
      betaM*((x[2])/(sum(x[-c(4,5)])))
    }
  } else {
    colnames(x) <- c("Sh", "Ih", "Rh", "Sm", "Im", paste("T", 1:n, sep=""))
    gam <- function(x){
      betaM*((x[2]+ psi[1:n] %*% x[6:(5+n)])/(sum(x[-c(4,5)])))
    }
  }
  change <- matrix(c(LamH, 0, 0, LamM, 0, 
                     0, 0, 0, 0, 0, 
                     0, 1/step -MuH-delS-sigma, sigma, 0, 0, 
                     omega, 0, 1/step - MuH - omega, 0, 0, 
                     0, 0, 0, 0, 0, 
                     0, 0, 0, 0, 1/step -MuM), nrow=6, byrow=T)
  if(n != 0){
    change2 <- matrix(0, nrow=n, ncol=n)
    for(i in 1:n-1){
      change2[i,i] <- 1/step - MuH - delS - alpha[[i+1]]
      change2[i, i+1] <- alpha[[i+1]]
    }
    change2[n,n] <- 1/step - MuH - delS - alpha[[n+1]]
    change <- rbind(cbind(change, matrix(0, nrow=6, ncol=n)), cbind(matrix(0, nrow=n, ncol=5), change2))
    change[3,2] <- 1/step -MuH-delS-sigma - alpha[[1]]
    change[3,6] <- alpha[[1]]
    change[6+n,1] <- alpha[[n]]
  }
  for(i in 2:iterations){
    lam = betaH*(x[i-1,5]/(sum(x[i-1,-c(4,5)]))) 
    gamm = gam(x[i-1,])
    change[2,1] <- 1/step - MuH - lam
    change[2,2] <- lam
    change[5,4] <- 1/step - MuM - gamm
    change[5,5] <- gamm
    x[i,] <- c(1,x[i-1,]) %*% (change*step)
    if(sum(abs((x[i-1,] - x[i,]))) < conv){
      print(paste("Convergence Successful in ", i, "iterations."))
      break
    }
  }
  rownames(x) <- seq(1, nrow(x))
  return(x[1:i,])
}

colH <- c("#08306B", "#08519C", "#2171B5", "#4D004B", "#810F7C", "#88419D", "#8C6BB1", "#8C96C6", "#67A9CF", "#3690C0", "#02818A", "#016C59", "#014636")
colM <- c("#67000D","#A50F15")
legH <- c("Susceptible Humans", "Infected Humans", "Recovered Humans", paste("Treatment", 1:10, sep=""))
legM <- c("Susceptible Mosquitoes","Infected Mosquitoes")


# Define server logic to plot various variables against mpg ----
server <- function(input, output) {

  time <- reactive(seq(0, input$time1, by=input$step))
  net <- reactive(c(rep(input$poph,3), rep(input$popm,2), rep(input$poph,input$bins)))
  results2N <- reactive(
    fthreeN(first=as.numeric(unlist(strsplit(input$vec3,",")))[1:(input$bins+5)]*net(), 
            iterations=length(time()), n=input$bins, conv=1*exp(-16),
            LamH=input$LamH, LamM=input$LamM, betaM=input$betaM, betaH=input$betaH,
            step=input$step, MuH=input$MuH, delS=input$delS, sigma=input$sigma,
            omega=input$omega, MuM=input$MuM, poph=input$poph, popm=input$popm,
            alpha=as.numeric(unlist(strsplit(input$vec1,",")))[1:(input$bins+1)], 
            psi=as.numeric(unlist(strsplit(input$vec2,",")))[1:input$bins])
  )
  
  output$HumanPlot <- renderPlot({
    plot(x=time(), results2N()[,1], type="l", lwd=2, col="blue", xlab="Time", ylab="Humans", ylim=c(0,max(results2N()[,-c(4:5)])))
    lines(x=time(), results2N()[,2], col=colH[2])
    lines(x=time(), results2N()[,3], col=colH[3])
    if(input$bins > 0){
      for(i in 6:(input$bins+5)){
        lines(x=time(), results2N()[,i], col=colH[(i-2)])
      }
    }
    legend(225,0.9*max(results2N()[,-c(4:5)]), legH[1:(3+input$bins)], lty=1, col=colH[1:(3+input$bins)], lwd=c(2, rep(1, (input$bins+2))))
  })
  
  # Compute the formula text ----
  # This is in a reactive expression since it is shared by the
  # output$caption and output$mpgPlot functions
  formulaText <- reactive({
    paste("Human population with", input$bins, "treatment bins.")
     })
  
  # Return the formula text for printing as a caption ----
  output$caption1 <- renderText({
    formulaText()
  })


  output$caption2 <- renderText({
    "Mosquito Population"
  })
  
  
  # Generate a plot
  output$MosquitoPlot <- renderPlot({
    plot(x=time(), results2N()[,4], type="l", lwd=2, col=colM[1], xlab="Time", ylab="Mosquitos", ylim=c(0,max(results2N()[,4:5])))
    lines(x=time(), results2N()[,5], col=colM[2])
    legend(200, 0.9*max(results2N()[,4:5]), legM, lty=1, col=colM, lwd=c(2,1))
  })
  
  
  
}


#runApp("C:/Users/michael.rocha/OneDrive - West Point/desktop/shinyapp/app.R")

#alternatively, from within this document,
shinyApp(ui, server)

