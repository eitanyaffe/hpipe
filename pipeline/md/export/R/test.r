
library(shiny)

# Define UI for application that draws a histogram
ui = shinyUI(fluidPage(

  # Application title
  titlePanel("Hello Shiny!"),

  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("bins",
                  "Number of bins:",
                  min = 1,
                  max = 50,
                  value = 30),
        verbatimTextOutput("coords")
    ),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot",
                 click = clickOpts(id="plot_click"),
                 hover = hoverOpts(
                     id="plot_hover", delay=50, delayType=c("throttle"), clip=TRUE, nullOutside=TRUE),
                 brush = brushOpts(
                     id = "plot_brush", fill = "#9cf", stroke = "#036", opacity = 0.25,
                     delay = 300, delayType = c("debounce", "throttle"),
                     clip = TRUE, direction = c("xy", "x", "y"), resetOnNew = FALSE))
    )
  )
))

# Define server logic required to draw a histogram
server = shinyServer(function(input, output) {
  state = reactiveValues(x=0, y=0)

  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should be automatically
  #     re-executed when inputs change
  #  2) Its output type is a plot

  output$distPlot <- renderPlot({
      plot(1:10)
      points(state$x, state$y)
      grid()
  })

  output$coords <- renderPrint({
    cat(sprintf("x:%f y:%f\n", input$plot_hover$x, input$plot_hover$y))
  })

  observeEvent(input$plot_click, {
      state$x = input$plot_click$x
      state$y = input$plot_click$y
  })
})

rl=function() {
    source("test.r")
    shinyApp(server=server, ui=ui)
}
