####script for Salvador et al. Morphology analysis to make an interactive UMAP shiny app where you can see the shapes the UMAP dots correspond to
####created by Madeline Melzer 20230915, last update by Madeline Melzer on 20240806

library(tidyverse)
library(ggplot2)
library(shiny)
library(plotly)
library(htmlwidgets)
library(base64enc)
library(shinyjs)
library(viridis)
library(bslib)

set.seed(23)

##defining the directories
getwd()

dataDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/results/umapTables/"
plotDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/plots/shiny/"
imageURIDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/individuals/"


data = as_tibble(read.table(paste0(dataDirectory, "20240806_umap_allMetrics.csv"), stringsAsFactors=F, header = T, sep = ","))
data$age <- as.character(data$age)
data$animal <- as.character(data$animal)

data$file <- trimws(data$file) # Remove any leading or trailing whitespace


continuous_vars <- c('area',
                     'perimeter',
                     'major_axis_length',
                     'minor_axis_length',
                     'eccentricity',
                     'solidity',
                     'equivalent_diameter_area',
                     'area_bbox',
                     'area_convex',
                     'extent',
                     'feret_diameter_max',
                     'perimeter_crofton',
                     'moments_hu.0',
                     'moments_hu.1',
                     'moments_hu.2',
                     'moments_hu.3',
                     'moments_hu.4',
                     'moments_hu.5',
                     'moments_hu.6',
                     'circularity',
                     'aspectRatio',
                     'flattening',
                     'averagePointRadius',
                     'pointRadiusExtremeDifference',
                     'pointRadiusRMSDifference',
                     'areaEquivalentEllipse',
                     'perimeterEquivalentEllipse',
                     'dArea',
                     'dPerimeter',
                     'dAreaConvex',
                     'xInterceptCount',
                     'concavityCount')
discrete_vars <- c("file", "label", "age", "sex", "animal", "location")

variable_names <- c("file", "label", "age", "sex", "animal", "location",
                    'area',
                    'perimeter',
                    'major_axis_length',
                    'minor_axis_length',
                    'eccentricity',
                    'solidity',
                    'equivalent_diameter_area',
                    'area_bbox',
                    'area_convex',
                    'extent',
                    'feret_diameter_max',
                    'perimeter_crofton',
                    'moments_hu.0',
                    'moments_hu.1',
                    'moments_hu.2',
                    'moments_hu.3',
                    'moments_hu.4',
                    'moments_hu.5',
                    'moments_hu.6',
                    'circularity',
                    'aspectRatio',
                    'flattening',
                    'averagePointRadius',
                    'pointRadiusExtremeDifference',
                    'pointRadiusRMSDifference',
                    'areaEquivalentEllipse',
                    'perimeterEquivalentEllipse',
                    'dArea',
                    'dPerimeter',
                    'dAreaConvex',
                    'xInterceptCount',
                    'concavityCount')

# Convert image to Data URI function
image_to_dataURI <- function(image_path) {
  # Read the binary image data
  img_data <- readBin(image_path, what = raw(), n = file.info(image_path)$size)
  
  # Convert to base64
  base64_img <- base64encode(img_data)
  
  # Construct the Data URI without the img tag
  data_uri <- paste0("data:image/png;base64,", base64_img)
  return(data_uri)
}

# Use the function to generate Data URIs for your images
data$image_dataURI <- mapply(function(f, l) {
  image_path <- paste0(imageURIDirectory, sub(".tif$", "", basename(f)), "_", l, ".png")
  image_to_dataURI(image_path)
}, data$file, data$label, SIMPLIFY = TRUE)




server <- function(input, output) {

  data$custom_info <- paste("<br>", "Age: ", data$age, " wk", "<br>",
                            "Location: ", data$location, "<br>",
                            "Sex: ", data$sex, "<br>",
                            "Animal: ", data$animal, "<br>")

  # Reactive variable to store the URI of the clicked point
  clickedURI <- reactiveVal(NULL)

  output$umapPlot <- renderPlotly({
    
    p <- ggplot(data, aes_string(x = "UMAP1", y = "UMAP2", color = input$color_var, customdata = "custom_info", ids = "image_dataURI"))
    p <- p +
      geom_point(size = 3) +
      theme_classic(base_size = 18) +
      coord_fixed(ratio = 1) +
      labs(x = "UMAP1", y = "UMAP2")
    
    # Conditional check for scale
    if (input$color_var %in% continuous_vars) {
      p <- p + scale_color_viridis_c()
    } else {
      p <- p + scale_color_viridis_d()
    }
    
    p_plotly = ggplotly(p, tooltip = "customdata")

    p_plotly <- onRender(p_plotly, "
       function(el) {
          el.on('plotly_hover', function(data) {
            // Adjust hover handling here
            // This can include more precise determination of the hovered point
          });
       
          el.on('plotly_click', function(data) {
             var uri = el.data[0].ids[data.points[0].pointNumber];
             Shiny.setInputValue('clicked_image_uri', uri); // Send the URI to Shiny

             var customInfo = data.points[0].customdata;
             $('.hoverlayer .hovertext').html(customInfo);
          });
       }
    ")

    return(p_plotly)
  })

  # Display the clicked image within the app
  output$clickedImage <- renderUI({
    current_data <- clicked_image_data()
    if (!is.null(current_data$uri)) {
      # Exclude the 'uri' from the displayed data
      displayed_data <- current_data[!names(current_data) %in% "uri"]
      
      # Generate the UI elements
      ui_elements <- lapply(names(displayed_data), function(name) {
        tags$p(paste(name, ":", displayed_data[[name]]))
      })
      
      # Create a layout with image in one column and metadata in another column
      fluidRow(
        column(6, tags$img(src = current_data$uri, width = "100%")), # Image column
        column(6, ui_elements)  # Metadata column
      )
    }
  })
  
  # Reactive values to store the image URI and associated metadata
  clicked_image_data <- reactiveVal(list(uri=NULL, age=NULL, sex=NULL, location=NULL))
  
  # On image click, update the reactive values with the current data
  observeEvent(input$clicked_image_uri, {
    image_index <- which(data$image_dataURI == input$clicked_image_uri)
    # Create a named list dynamically
    clicked_values <- sapply(variable_names, function(var) {
      data[[var]][image_index]
    }, USE.NAMES = TRUE)
    # Add the image URI
    clicked_values$uri <- input$clicked_image_uri
    # Store the data
    clicked_image_data(clicked_values)
  })

}


ui <- fluidPage(
  tags$head(
    # Add custom CSS for background color and font
    tags$style(HTML("
      body {
        background-color: #FFFFFF;
        font-family: 'Arial', sans-serif;
        font-weight: bold;
      }
      .shiny-title-panel {
        text-align: center;
      }
      h1, h2 {
        text-align: center;
      }
    "))
  ),
  
  titlePanel("UMAP of nuclear features"),
  
  #Subheader
  div(
    style = "font-size: 24px; text-align:center;",   # Adjust this value for desired font size
    "Salvador, J., … Melzer, M.E.,... Goyal, Y., … Iruela Arispe, M.L., 2023", 
    tags$br()
    
  ), 
  
  div(
    style = "font-size: 18px; text-align:center;",
    "Created by MEM on 20230915", 
    tags$br(), 
    "Last update by MEM on 20231117"
  ),
  
  selectInput("color_var",
              "Choose a variable to color by:",
              choices = c('age', 
                          'sex', 
                          'animal', 
                          'location', 
                          'area',
                          'perimeter',
                          'major_axis_length',
                          'minor_axis_length',
                          'eccentricity',
                          'solidity',
                          'equivalent_diameter_area',
                          'area_bbox',
                          'area_convex',
                          'extent',
                          'feret_diameter_max',
                          'perimeter_crofton',
                          'moments_hu.0',
                          'moments_hu.1',
                          'moments_hu.2',
                          'moments_hu.3',
                          'moments_hu.4',
                          'moments_hu.5',
                          'moments_hu.6',
                          'circularity',
                          'aspectRatio',
                          'flattening',
                          'averagePointRadius',
                          'pointRadiusExtremeDifference',
                          'pointRadiusRMSDifference',
                          'areaEquivalentEllipse',
                          'perimeterEquivalentEllipse',
                          'dArea',
                          'dPerimeter',
                          'dAreaConvex',
                          'xInterceptCount',
                          'concavityCount')),
  
  # Another Subheader
  #tags$h2("Visualization and Details:"),

  
  fluidRow(
    column(6, plotlyOutput("umapPlot", width = "1000px", height = "850px")), # This takes up 8/12 of the row's width
    column(5, uiOutput("clickedImage"))   # This takes up 4/12 of the row's width
  )
)

shinyApp(ui = ui, server = server)




#########Revised


server <- function(input, output) {
  
  data$custom_info <- paste("<br>", "Age: ", data$age, " wk", "<br>",
                            "Location: ", data$location, "<br>",
                            "Sex: ", data$sex, "<br>",
                            "Animal: ", data$animal, "<br>")
  
  # Reactive variable to store the URI of the clicked point
  clickedURI <- reactiveVal(NULL)
  
  output$umapPlot <- renderPlotly({
    # ... [plotting code remains the same]
    
    # Modify the onRender script
    p_plotly <- onRender(p_plotly, "
       function(el) {
          el.on('plotly_hover', function(data) {
            // Adjust hover handling here
            // This can include more precise determination of the hovered point
          });
          el.on('plotly_click', function(data) {
             // Adjust click handling here
             // Ensure that the clicked point's data is accurately captured
          });
       }
    ")
    
    
    
    return(p_plotly)
  })
  
  # ... [rest of your code remains the same]
}






