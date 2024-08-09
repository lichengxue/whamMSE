library(shiny)
library(DiagrammeR)

# Define global variable
mean_vals <- NULL

# Define UI for the app
ui <- fluidPage(
  titlePanel("Specify Movement Rates for Fish"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("n_stocks", "Number of Stocks:", min = 1, max = 10, value = 2),
      numericInput("n_seasons", "Number of Seasons:", min = 1, max = 4, value = 4),
      numericInput("n_regions", "Number of Regions:", min = 2, max = 10, value = 2),
      
      uiOutput("movementInputs"),
      
      actionButton("generate", "Generate Movement Matrix"),
      downloadButton("downloadData", "Download Movement Matrix")
    ),
    
    mainPanel(
      textOutput("instructions"),
      textOutput("warning"),
      grVizOutput("movementDiagram"), # Move the diagram here
      verbatimTextOutput("meanValsArray")
    )
  )
)

# Define server logic for the app
server <- function(input, output, session) {
  observe({
    # Generate dynamic UI inputs for movement rates
    n_regions <- input$n_regions
    movement_inputs <- list()
    
    for (r in 1:n_regions) {
      for (k in 1:(n_regions - 1)) {
        rr <- if (k < r) k else k + 1
        movement_inputs <- append(movement_inputs, list(
          numericInput(paste0("move_", r, "_", k),
                       paste0("Movement rate from region ", r, " to region ", rr, ":"),
                       value = 0, min = 0, max = 1, step = 0.01)
        ))
      }
    }
    
    output$movementInputs <- renderUI({
      do.call(tagList, movement_inputs)
    })
  })
  
  output$instructions <- renderText({
    n_stocks <- input$n_stocks
    n_regions <- input$n_regions
    if (is.null(n_stocks) || is.null(n_regions)) {
      return("")  # Return an empty string to avoid evaluating missing values
    }
    required_values <- n_stocks * (n_regions * (n_regions - 1))
    paste("You need to specify", required_values, "movement values for", n_stocks, "stocks and", n_regions, "regions.")
  })
  
  output$warning <- renderText({
    n_stocks <- input$n_stocks
    n_regions <- input$n_regions
    if (is.null(n_stocks) || is.null(n_regions)) {
      return("")  # Avoid rendering warning if inputs are not yet initialized
    }
    if (n_stocks != n_regions) {
      paste("Warning: The number of stocks (", n_stocks, ") is not equal to the number of regions (", n_regions, ").")
    } else {
      ""
    }
  })
  
  mean_vals <- reactiveVal(NULL)
  
  observeEvent(input$generate, {
    # Generate the mean_vals array
    n_stocks <- input$n_stocks
    n_seasons <- input$n_seasons
    n_regions <- input$n_regions
    
    mean_vals_array <- array(0, dim = c(n_stocks, n_seasons, n_regions, n_regions - 1))
    
    for (r in 1:n_regions) {
      k_index <- 1
      for (rr in 1:n_regions) {
        if (r != rr) {
          input_value <- as.numeric(input[[paste0("move_", r, "_", k_index)]])
          if (is.na(input_value)) {
            input_value <- 0 # Set a default value in case input is NA
          }
          mean_vals_array[,,r,k_index] <- input_value
          k_index <- k_index + 1
        }
      }
    }
    
    mean_vals(mean_vals_array)
    
    output$meanValsArray <- renderPrint({ mean_vals_array })
    
    # Create the movement diagram using DiagrammeR
    diagram <- create_graph() %>%
      add_global_graph_attrs(attr = "layout", value = "dot", attr_type = "graph") %>%
      add_global_graph_attrs(attr = "rankdir", value = "LR", attr_type = "graph")
    
    for (r in 1:n_regions) {
      diagram <- diagram %>% add_node(label = paste("Region", r))
    }
    
    for (r in 1:n_regions) {
      k_index <- 1
      for (rr in 1:n_regions) {
        if (r != rr) {
          input_value <- as.numeric(input[[paste0("move_", r, "_", k_index)]])
          if (input_value > 0) {  # Only add arrows where movement is greater than 0
            diagram <- diagram %>%
              add_edge(from = r, to = rr, edge_aes = edge_aes(label = paste0(input_value)))
          }
          k_index <- k_index + 1
        }
      }
    }
    
    output$movementDiagram <- renderGrViz({
      grViz(DiagrammeR::generate_dot(diagram))
    })
    
    # Save to global environment
    isolate({
      assign("mean_vals", mean_vals_array, envir = .GlobalEnv)
    })
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("mean_vals", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      mean_vals_array <- mean_vals()
      mean_vals_flattened <- as.data.frame.table(mean_vals_array)
      colnames(mean_vals_flattened) <- c("Stock", "Season", "From_Region", "To_Region_Index", "Value")
      write.csv(mean_vals_flattened, file, row.names = FALSE)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)