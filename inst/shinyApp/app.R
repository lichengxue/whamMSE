library(shiny)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(shinyjs)

# Define UI for the app
ui <- fluidPage(
  useShinyjs(),  # Enable shinyjs
  
  titlePanel("Specify Movement Rates for Fish"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("n_stocks", "Number of Stocks:", min = 1, max = 100, value = 2),
      numericInput("n_seasons", "Number of Seasons:", min = 1, max = 100, value = 4),
      numericInput("n_regions", "Number of Regions:", min = 2, max = 100, value = 2),
      numericInput("fracyr_spawn", "Fraction of Year for Spawning:", min = 0, max = 1, value = 0.5, step = 0.1),
      
      checkboxGroupInput("canMoveCheckbox", "Can Move in Seasons:", choices = c(1, 2, 3, 4), selected = c(1, 2, 3, 4)),
      checkboxGroupInput("mustMoveCheckbox", "Must Move in Seasons:", choices = c(1, 2, 3, 4), selected = 1),
      
      uiOutput("moveEffectOptions"),
      uiOutput("moveSigmaInput"),
      uiOutput("moveRhoAInput"),
      uiOutput("moveRhoYInput"),
      checkboxInput("use_prior", "Use Prior?", value = FALSE),
      uiOutput("priorSigmaInput"),
      uiOutput("movementInputs"),
      
      actionButton("generate", "Generate Movement Matrix"),
      downloadButton("downloadMovementMatrix", "Download Movement Matrix"),
      downloadButton("downloadMovementDiagram", "Download Movement Diagram"),
      actionButton("demoButton", "Demo"),
      actionButton("restartButton", "Restart"),
      actionButton("exitButton", "Exit")
    ),
    
    mainPanel(
      uiOutput("warning"),
      grVizOutput("movementDiagram"),
      plotOutput("seasonDiagram"),
      verbatimTextOutput("outputList")
    )
  )
)

# Define server logic for the app
server <- function(input, output, session) {
  
  diagram_reactive <- reactiveVal(NULL)  # Create a reactive value to store the diagram
  output_list_reactive <- reactiveVal(NULL)  # Create a reactive value to store the output_list
  
  observe({
    n_seasons <- input$n_seasons
    fracyr_spawn <- input$fracyr_spawn
    
    # Calculate spawning season
    spawn_season <- ceiling(fracyr_spawn * n_seasons)
    
    # Update Can Move and Must Move checkboxes based on spawning season
    updateCheckboxGroupInput(session, "canMoveCheckbox", choices = 1:n_seasons, selected = setdiff(1:n_seasons, spawn_season))
    updateCheckboxGroupInput(session, "mustMoveCheckbox", choices = 1:n_seasons, selected = spawn_season)
    
    # Generate movement effect options
    output$moveEffectOptions <- renderUI({
      selectInput("movementEffect", "Select Movement Random Effect",
                  choices = list("Constant" = "constant", "IID (Age)" = "iid_a",
                                 "IID (Year)" = "iid_y", "AR1 (Age)" = "ar1_a",
                                 "AR1 (Year)" = "ar1_y"), selected = "constant")
    })
    
    # Conditionally show move.sigma input
    output$moveSigmaInput <- renderUI({
      if (!is.null(input$movementEffect) && input$movementEffect %in% c("iid_a", "iid_y", "ar1_a", "ar1_y")) {
        numericInput("moveSigma", "Movement Sigma:", value = 0.5, step = 0.1)
      }
    })
    
    # Conditionally show move.rho.a input
    output$moveRhoAInput <- renderUI({
      if (!is.null(input$movementEffect) && input$movementEffect == "ar1_a") {
        numericInput("moveRhoA", "Movement Rho (Age):", value = 0.5, step = 0.1)
      }
    })
    
    # Conditionally show move.rho.y input
    output$moveRhoYInput <- renderUI({
      if (!is.null(input$movementEffect) && input$movementEffect == "ar1_y") {
        numericInput("moveRhoY", "Movement Rho (Year):", value = 0.5, step = 0.1)
      }
    })
    
    # Conditionally show prior.sigma input
    output$priorSigmaInput <- renderUI({
      if (input$use_prior) {
        numericInput("priorSigma", "Prior Sigma:", value = 0.2, step = 0.1)
      }
    })
    
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
    
    # Generate season diagram
    output$seasonDiagram <- renderPlot({
      plot(1:n_seasons, rep(1, n_seasons), xlim = c(0, n_seasons + 1), ylim = c(0, 2),
           xaxt = "n", yaxt = "n", xlab = "", ylab = "", type = "n")
      rect(1:n_seasons - 0.4, 0.5, 1:n_seasons + 0.4, 1.5,
           col = ifelse(1:n_seasons == spawn_season, "#35B779FF", "#FDE725FF"))
      text(1:n_seasons, rep(1, n_seasons), labels = ifelse(1:n_seasons == spawn_season, "Spawning\n Season", "Offspawning\n Season"), cex = 0.8)
      axis(1, at = 1:n_seasons, labels = paste("Season", 1:n_seasons))
    })
    
    mean_vals <- reactiveVal(NULL)
    
    observeEvent(input$generate, {
      # Lock all inputs except Restart and Exit
      lapply(names(input), function(id) {
        if (id != "restartButton" && id != "exitButton") {
          shinyjs::disable(id)
        }
      })
      
      # Generate the mean_vals array and other outputs
      n_stocks <- input$n_stocks
      n_seasons <- input$n_seasons
      n_regions <- input$n_regions
      use_prior <- input$use_prior
      prior_sigma <- if (use_prior) input$priorSigma else NULL
      
      # Assuming home_region for each stock is stock index (e.g., Stock 1 home in Region 1, Stock 2 home in Region 2, etc.)
      home_region <- 1:n_stocks
      
      mean_vals_array <- array(0, dim = c(n_stocks, n_seasons, n_regions, n_regions - 1))
      
      error_detected <- FALSE  # Initialize error flag
      
      for (r in 1:n_regions) {
        k_index <- 1
        sum_of_movements <- 0  # Initialize the sum for each region
        
        for (rr in 1:n_regions) {
          if (r != rr) {
            input_id <- paste0("move_", r, "_", k_index)
            input_value <- as.numeric(input[[input_id]])
            
            if (is.na(input_value)) {
              input_value <- 0  # Default to zero if input is missing
            }
            
            # Assign input_value to mean_vals_array
            mean_vals_array[, , r, k_index] <- input_value
            
            sum_of_movements <- sum_of_movements + input_value
            k_index <- k_index + 1
          }
        }
        
        # Display an error if the sum exceeds 1
        if (sum_of_movements > 1) {
          output$warning <- renderUI({
            HTML(paste0("<span style='color: red; font-size: 20px;'>Warning: The sum of movement rates for Region ", r, " exceeds 1. Please adjust the values.</span>"))
          })
        }
      }
      
      mean_vals(mean_vals_array)
      
      # Example calculation for other outputs based on selections
      can_move_array <- array(1, dim = c(n_stocks, n_seasons, n_regions, n_regions))
      must_move_array <- array(0, dim = c(n_stocks, n_seasons, n_regions))
      
      must_move_seasons <- as.integer(input$mustMoveCheckbox)
      
      for (stock in 1:n_stocks) {
        for (season in 1:n_seasons) {
          for (r in 1:n_regions) {
            for (rr in 1:n_regions) {
              if (!(season %in% input$canMoveCheckbox)) {
                can_move_array[stock, season, r, rr] <- 0
              }
            }
            if (season %in% input$mustMoveCheckbox && r != home_region[stock]) {
              must_move_array[stock, season, r] <- 1
            }
          }
        }
      }
      
      mean_model_matrix <- matrix("constant", nrow = n_regions, ncol = n_regions - 1)
      
      year_re_matrix <- matrix("none", nrow = n_regions, ncol = n_regions - 1)
      age_re_matrix <- matrix("none", nrow = n_regions, ncol = n_regions - 1)
      
      cor_vals <- array(0, dim = c(n_stocks, n_seasons, n_regions, n_regions - 1, 2))  # Corrected dimensions
      if (input$movementEffect == "ar1_a") {
        cor_vals[, , , , 1] <- input$moveRhoA  # Assign age correlation
        age_re_matrix <- matrix("ar1", nrow = n_regions, ncol = n_regions - 1)
      }
      if (input$movementEffect == "ar1_y") {
        cor_vals[, , , , 2] <- input$moveRhoY  # Assign year correlation
        year_re_matrix <- matrix("ar1", nrow = n_regions, ncol = n_regions - 1)
      }
      
      sigma_vals <- if (input$movementEffect %in% c("iid_a", "iid_y", "ar1_a", "ar1_y")) input$moveSigma else NULL
      
      use_prior_array <- array(0, dim = c(n_stocks, n_seasons, n_regions, n_regions - 1))
      if (use_prior) {
        use_prior_array[, 1, , ] <- 1  # Set the first season to 1, others remain 0
      }
      
      output_list <- list(stock_move = rep(TRUE, n_regions),
                          separable = TRUE,
                          must_move = must_move_array,
                          can_move = can_move_array,
                          mean_vals = mean_vals_array,
                          mean_model = mean_model_matrix,
                          year_re = year_re_matrix,
                          age_re = age_re_matrix,
                          cor_vals = cor_vals,
                          sigma_vals = sigma_vals,
                          use_prior = use_prior_array,
                          prior_sigma = if (!is.null(prior_sigma)) array(prior_sigma, dim = c(n_stocks, n_seasons, n_regions, n_regions - 1)) else NULL)
      
      output_list_reactive(output_list)  # Store output_list in the reactive value
      
      output$outputList <- renderPrint({ output_list })
      
      # Create the movement diagram using DiagrammeR
      diagram <- create_graph() %>%
        add_global_graph_attrs(attr = "layout", value = "dot", attr_type = "graph") %>%
        add_global_graph_attrs(attr = "rankdir", value = "LR", attr_type = "graph")
      
      for (r in 1:n_regions) {
        diagram <- diagram %>% add_node(label = paste("Region\n", r))
      }
      
      for (r in 1:n_regions) {
        k_index <- 1
        for (rr in 1:n_regions) {
          if (r != rr) {
            input_value <- as.numeric(input[[paste0("move_", r, "_", k_index)]])
            if (input_value > 0) {
              diagram <- diagram %>%
                add_edge(from = r, to = rr, edge_aes = edge_aes(label = paste0(input_value)))
            }
            k_index <- k_index + 1
          }
        }
      }
      
      diagram_reactive(diagram)  # Store the diagram in the reactive value
      
      output$movementDiagram <- renderGrViz({
        grViz(DiagrammeR::generate_dot(diagram_reactive()))
      })
      
      # Save to global environment
      isolate({
        assign("output_list", output_list, envir = .GlobalEnv)
      })
    })
    
    # Demo Button: Pre-set values for demo
    observeEvent(input$demoButton, {
      updateNumericInput(session, "n_stocks", value = 4)
      updateNumericInput(session, "n_regions", value = 4)
      updateNumericInput(session, "n_seasons", value = 4)
      updateNumericInput(session, "fracyr_spawn", value = 0.5)
      
      for (r in 1:4) {
        k_index <- 1
        remaining_sum <- 1
        for (rr in 1:4) {
          if (r != rr) {
            random_value <- round(runif(1, 0, remaining_sum), 2)
            updateNumericInput(session, paste0("move_", r, "_", k_index), value = random_value)
            remaining_sum <- remaining_sum - random_value
            k_index <- k_index + 1
          }
        }
      }
    })
    
    # Restart Button: Reload the app
    observeEvent(input$restartButton, {
      session$reload()
    })
    
    # Exit the app
    observeEvent(input$exitButton, {
      stopApp()
    })
  })
  
  output$downloadMovementDiagram <- downloadHandler(
    filename = function() {
      paste("movement_diagram.png", sep = "")
    },
    content = function(file) {
      export_svg(grViz(DiagrammeR::generate_dot(diagram_reactive()))) %>%
        charToRaw %>%
        rsvg_png(file)
    }
  )
  
  output$downloadMovementMatrix <- downloadHandler(
    filename = function() {
      paste("movement_matrix.rds", sep = "")
    },
    content = function(file) {
      saveRDS(output_list_reactive(), file)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)