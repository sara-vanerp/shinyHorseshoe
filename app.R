##### Shiny Regularized Horseshoe #####
## Author: Sara van Erp 

library(shiny)
library(LaplacesDemon)
library(ggplot2)

# Define UI 
ui <- fluidPage(

    # Application title
    titlePanel(HTML(paste(h1("Regularized horseshoe visualization"), h5("Developed by:", a("Sara van Erp", href = "https://saravanerp.com") )))),
    br(),
    
    fluidRow(
      
      column(3,
             htmlOutput("general"),
             br(),
             numericInput("ndraws",
                          "Number of draws from the prior",
                          value = 5000),
             numericInput("xlow",
                          "Lower limit x-axis",
                          value = -5),
             numericInput("xhigh",
                          "Upper limit x-axis",
                          value = 5)
             ),
      
      column(3,
             htmlOutput("prior1"),
             br(),
             numericInput("nu1_1",
                          HTML("&nu;<sub>1</sub>"),
                          value = 1),
             numericInput("nu2_1",
                          HTML("&nu;<sub>2</sub>"),
                          value = 1),
             numericInput("nu3_1",
                          HTML("&nu;<sub>3</sub>"),
                          value = 4),
             numericInput("s_1", 
                          HTML("<em>s</em>"),
                          value = 1),
             numericInput("lambda20_1", 
                          HTML("&lambda;<sup>2</sup><sub>0</sub>"),
                          value = 1)
             ),
      
      column(3,
             htmlOutput("prior2"),
             br(),
             numericInput("nu1_2",
                          HTML("&nu;<sub>1</sub>"),
                          value = 1),
             numericInput("nu2_2",
                          HTML("&nu;<sub>2</sub>"),
                          value = 1),
             numericInput("nu3_2",
                          HTML("&nu;<sub>3</sub>"),
                          value = 4),
             numericInput("s_2", 
                          HTML("<em>s</em>"),
                          value = 1),
             numericInput("lambda20_2", 
                          HTML("&lambda;<sup>2</sup><sub>0</sub>"),
                          value = 1)
      ),
      
      column(3,
             htmlOutput("prior3"),
             br(),
             numericInput("nu1_3",
                          HTML("&nu;<sub>1</sub>"),
                          value = 1),
             numericInput("nu2_3",
                          HTML("&nu;<sub>2</sub>"),
                          value = 1),
             numericInput("nu3_3",
                          HTML("&nu;<sub>3</sub>"),
                          value = 4),
             numericInput("s_3", 
                          HTML("<em>s</em>"),
                          value = 1),
             numericInput("lambda20_3", 
                          HTML("&lambda;<sup>2</sup><sub>0</sub>"),
                          value = 1)
      )
    ),
    
    hr(),

    mainPanel(
      
    column(12, 
           plotOutput("priorPlot")
           )
    )
)

# Define server logic 
server <- function(input, output) {

    output$priorPlot <- renderPlot({
      # function to sample from the prior
      priorfun <- function(ndraws, nu1, nu2, nu3, s, lambda20){
        reghs <- rep(NA, ndraws)
        for(i in 1:ndraws){
          c2 <- rinvgamma(1, shape=nu3/2, scale=nu3*s^2/2)
          lambda <- rhalft(1, scale=lambda20, nu = nu1)
          tau <- rhalft(1, scale=1, nu = nu2)
          tau2_tilde <- c2 * tau^2/(c2 + lambda^2*tau^2)
          reghs[i] <- rnorm(1, 0, sqrt(lambda^2*tau2_tilde))
        }
        return(reghs)
      }
      
      reghs1 <- priorfun(ndraws = input$ndraws, 
                         nu1 = input$nu1_1,
                         nu2 = input$nu2_1,
                         nu3 = input$nu3_1,
                         s = input$s_1,
                         lambda20 = input$lambda20_1)
      df1 <- data.frame("draws" = reghs1, "Prior" = "Prior setting 1")
      reghs2 <- priorfun(ndraws = input$ndraws, 
                         nu1 = input$nu1_2,
                         nu2 = input$nu2_2,
                         nu3 = input$nu3_2,
                         s = input$s_2,
                         lambda20 = input$lambda20_2)
      df2 <- data.frame("draws" = reghs2, "Prior" = "Prior setting 2")
      reghs3 <- priorfun(ndraws = input$ndraws, 
                         nu1 = input$nu1_3,
                         nu2 = input$nu2_3,
                         nu3 = input$nu3_3,
                         s = input$s_3,
                         lambda20 = input$lambda20_3)
      df3 <- data.frame("draws" = reghs3, "Prior" = "Prior setting 3")
      
      df <- rbind.data.frame(df1, df2, df3)
      
      # plot
      ggplot(df, aes(x = draws, group = Prior)) +
        geom_density(aes(colour = Prior)) +
        coord_cartesian(xlim = c(input$xlow, input$xhigh)) +
        theme_bw(base_size = 14) + theme(legend.title = element_blank()) +
        xlab("") + ylab("Density") + ggtitle("Comparison of various regularized horseshoe specifications")
    })
    
    output$general <- renderText("<font size =+2>General</font>")
    output$prior1 <- renderText("<font size =+2>Prior setting 1</font>")
    output$prior2 <- renderText("<font size =+2>Prior setting 2</font>")
    output$prior3 <- renderText("<font size =+2>Prior setting 3</font>")
}

# Run the application 
shinyApp(ui = ui, server = server)
