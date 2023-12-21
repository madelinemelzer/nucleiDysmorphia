####script for Melzer et al. Morphology analysis to make an interactive UMAP shiny app where you can see the shapes the points correspond to. 
####created by Madeline Melzer 20230915
#works when ran in RStudio but not DataSpell?

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

## PC:
#dataDirectory <- "R:Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/output_metrics/UMAPs/"
#plotDirectory <- "R:Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/plots/"

## Mac:
dataDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/output_metrics/UMAPs/"
plotDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/plots/"
imageURIDirectory = "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/nucleiSegmented/individuals/png/"


data = as_tibble(read.table(paste0(dataDirectory, "20231221_umap_allAges.csv"), stringsAsFactors=F, header = T, sep = ","))
data$age <- as.character(data$age)
data$animal <- as.character(data$animal)

data$file <- trimws(data$file) # Remove any leading or trailing whitespace


continuous_vars <- c('eccentricity', 'perimeter', 'dArea', 'area', 'major_axis_length', 'minor_axis_length', 'r_avg', 'rms_avg', 'rms_subsampled', 'convex_area', 'x_intercepts', 'concavity_ct', 'circularity', 'solidity', 'equivalent_diameter_area', 'area_bbox', 'area_convex', 'area_filled', 'extent', 'feret_diameter_max', 'minkowskiDimension', 'compactness', 'roundness', 'convexity', 'aspect_ratio', 'perimeter_crofton', 'concavity_ct.log')
discrete_vars <- c("age", "sex", "animal", "location")

variable_names <- c('age', 'sex', 'animal', 'location', 'eccentricity', 
                    'perimeter', 'dArea', 'area', 'major_axis_length', 
                    'minor_axis_length', 'r_avg', 'rms_avg', 'rms_subsampled', 
                    'convex_area', 'x_intercepts', 'concavity_ct', 'circularity', 
                    'solidity', 'equivalent_diameter_area', 'area_bbox', 
                    'area_convex', 'area_filled', 'extent', 'feret_diameter_max', 
                    'minkowskiDimension', 'compactness', 'roundness', 'convexity', 
                    'aspect_ratio', 'perimeter_crofton', 'concavity_ct.log')

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
    #if (input$color_var %in% c("age", "animal")) {
    #  p <- ggplot(data, aes_string(x = "UMAP1", y = "UMAP2", color = paste0("as.factor(", input$color_var, ")"), customdata = "custom_info", ids = "image_dataURI"))
    #} else {
    #  p <- ggplot(data, aes_string(x = "UMAP1", y = "UMAP2", color = input$color_var, customdata = "custom_info", ids = "image_dataURI"))
    #}
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
    "Last update by MEM on 20231023"
  ),
  
  selectInput("color_var",
              "Choose a variable to color by:",
              choices = c('age', 
                          'sex', 
                          'animal', 
                          'location', 
                          'eccentricity', 
                          'perimeter', 
                          'dArea', 
                          'area', 
                          'major_axis_length', 
                          'minor_axis_length', 
                          'r_avg', 
                          'rms_avg', 
                          'rms_subsampled', 
                          'convex_area', 
                          'x_intercepts', 
                          'concavity_ct', 
                          'circularity', 
                          'solidity', 
                          'equivalent_diameter_area', 
                          'area_bbox', 
                          'area_convex', 
                          'area_filled', 
                          'extent', 
                          'feret_diameter_max', 
                          'minkowskiDimension', 
                          'compactness', 
                          'roundness', 
                          'convexity', 
                          'aspect_ratio', 
                          'perimeter_crofton', 
                          'concavity_ct.log')),
  
  # Another Subheader
  #tags$h2("Visualization and Details:"),

  
  fluidRow(
    column(6, plotlyOutput("umapPlot", width = "1000px", height = "850px")), # This takes up 8/12 of the row's width
    column(5, uiOutput("clickedImage"))   # This takes up 4/12 of the row's width
  )
)



shinyApp(ui = ui, server = server)











































#OLD############################################################################################################################################
################################################################################################################################################
################################################################################################################################################




ui <- fluidPage(
  titlePanel("Interactive UMAP Plot"),
  plotlyOutput("umapPlot")
)

server <- function(input, output) {
  output$umapPlot <- renderPlotly({
    p <- ggplot(data, aes(x = UMAP1, y = UMAP2, color = location, text = paste("File Name:", file, "<br>Label:", label))) +
      geom_point(size = 3) +
      theme_classic(base_size = 18) +
      coord_fixed(ratio = 1)+
      labs(x = "UMAP1", y = "UMAP2")

    # Convert ggplot object to plotly object
    p_plotly = ggplotly(p, tooltip="text")
    p_plotly <- layout(p_plotly, autosize = F, width = 1000, height = 1000)

  })
}

shinyApp(ui = ui, server = server)


### With IMAGES of nuclei being pulled up
library(htmlwidgets)

dataDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/output_metrics/UMAPs/"
plotDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/plots/"

data = as_tibble(read.table(paste0(dataDirectory, "20230913_controlCombinedumapdf_relabeled_filteredAreas.csv"), stringsAsFactors=F, header = T, sep = ","))

data$file <- trimws(data$file) # Remove any leading or trailing whitespace

ui <- fluidPage(
  titlePanel("Interactive UMAP Plot"),
  plotlyOutput("umapPlot")
)

server <- function(input, output) {
  output$umapPlot <- renderPlotly({
    p <- ggplot(data, aes(x = UMAP1, y = UMAP2, color = location,
                          customdata = paste0("<img src='file:///Users/mem3579/dataspell/nucleidysmorphia/scripts/www/", sub(".png$", "", basename(file)), "_", label, ".png' width='200px'>",
                                              "Age: ", age, "<br>",
                                              "Location: ", location, "<br>",
                                              "Sex: ", sex, "<br>",
                                              "Animal: ", animal))) +
      geom_point(size = 3) +
      theme_classic(base_size = 18) +
      coord_fixed(ratio = 1) +
      labs(x = "UMAP1", y = "UMAP2")

    # Convert ggplot object to plotly object
    p_plotly = ggplotly(p, tooltip = "customdata")

    p_plotly <- onRender(p_plotly, "
       function(el) {
          el.on('plotly_hover', function(data) {
             var customInfo = data.points[0].customdata;
             $('.hoverlayer .hovertext').html(customInfo);
          });
       }
    ")

    return(p_plotly)

  })
}

shinyApp(ui = ui, server = server)


# Use the function to generate Data URIs for your images
data$image_dataURI <- mapply(function(f, l) {
  image_path <- paste0("/Users/mem3579/dataspell/nucleidysmorphia/scripts/www/", sub(".tif", "", basename(f)), "_", l, ".tif")
  image_to_dataURI(image_path)
}, data$file, data$label, SIMPLIFY = TRUE)

ui <- fluidPage(
  titlePanel("Interactive UMAP Plot"),
  plotlyOutput("umapPlot")
)

server <- function(input, output) {
  output$umapPlot <- renderPlotly({
    p <- ggplot(data, aes(x = UMAP1, y = UMAP2, color = location,
                          customdata = paste0("<img src='", image_dataURI, "' width='200px'><br>",
                                              "Age: ", age, "<br>",
                                              "Location: ", location, "<br>",
                                              "Sex: ", sex, "<br>",
                                              "Animal: ", animal))) +
      geom_point(size = 3) +
      theme_classic(base_size = 18) +
      coord_fixed(ratio = 1) +
      labs(x = "UMAP1", y = "UMAP2")

    # Convert ggplot object to plotly object
    p_plotly = ggplotly(p, tooltip = "customdata")

    p_plotly <- onRender(p_plotly, "
       function(el) {
          el.on('plotly_hover', function(data) {
             var customInfo = data.points[0].customdata;
             $('.hoverlayer .hovertext').html(customInfo);
          });
       }
    ")

    return(p_plotly)
  })
}


#tester using plotly_click

# Sample data for demonstration
ui <- fluidPage(
  titlePanel("Interactive UMAP Plot"),
  plotlyOutput("umapPlot")
)

server <- function(input, output, session) {

  observeEvent(event_data("plotly_click", "umapPlot"), {
    clicked_point <- event_data("plotly_click", "umapPlot")

    if(!is.null(clicked_point)) {
      point_index <- clicked_point$pointNumber + 1
      showModal(modalDialog(
        title = "Image",
        tags$img(src = data$image_dataURI[point_index], width = "100%")
      ))
    }
  })

  output$umapPlot <- renderPlotly({
    p <- plot_ly(data, x = ~UMAP1, y = ~UMAP2, color = ~location, ids = ~row.names(data),
                 type = "scatter", mode = "markers", marker = list(size = 10, opacity = 0.7)) %>%
      layout(hovermode = "closest", title = "Interactive UMAP Plot") %>%
      layout(coloraxis = list(colorscale = "Viridis"))

    # Register the plotly_click event
    event_register(p, "plotly_click")

    return(p)
  })

}

shinyApp(ui = ui, server = server)


#bypassing ggplot and using plot_ly directly

library(shiny)
library(plotly)
library(base64enc)
library(shiny)
library(plotly)
library(base64enc)

# Convert image to Data URI function
image_to_dataURI <- function(image_path) {
  # Read the binary image data
  img_data <- readBin(image_path, what = raw(), n = file.info(image_path)$size)

  # Convert to base64
  base64_img <- base64encode(img_data)

  # Return the Data URI
  data_uri <- paste0("data:image/png;base64,", base64_img)
  return(data_uri)
}

dataDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/output_metrics/UMAPs/"
plotDirectory <- "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/plots/"

#data = as_tibble(read.table(paste0(dataDirectory, "20230913_controlCombinedumapdf_relabeled_filteredAreas.csv"), stringsAsFactors=F, header = T, sep = ","))
#data$file <- trimws(data$file) # Remove any leading or trailing whitespace

# Use the function to generate Data URIs for your images
data$image_dataURI <- mapply(function(f, l) {
  image_path <- paste0("/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/nucleiSegmented/individuals/png/", sub(".tif$", "", basename(f)), "_", l, ".png")
  image_to_dataURI(image_path)
}, data$file, data$label, SIMPLIFY = TRUE)

ui <- fluidPage(
  titlePanel("Interactive UMAP Plot"),
  plotlyOutput("umapPlot")
)

server <- function(input, output) {
  output$umapPlot <- renderPlotly({
    plot_ly(data, x = ~UMAP1, y = ~UMAP2, color = ~location, ids = ~file, type = "scatter", mode = "markers",
            marker = list(size = 10, opacity = 1),
            hoverinfo = "text",
            hovertemplate = paste0("<img src='%{text}' width='200px' /><extra></extra>"),
            text = ~image_dataURI) %>%
      layout(hovermode = "closest", title = "Interactive UMAP Plot", coloraxis = list(colorscale = "Viridis"))
  })
}

shinyApp(ui = ui, server = server)


##### trying again with www


data$image_url <- paste0("/", sub(".tif$", ".png", basename(data$file)), "_", data$label)



ui <- fluidPage(
  titlePanel("Interactive UMAP Plot"),
  plotlyOutput("umapPlot")
)

server <- function(input, output) {
  output$umapPlot <- renderPlotly({
    plot_ly(data, x = ~UMAP1, y = ~UMAP2, color = ~location, text = ~image_url, ids = ~file,
            type = "scatter", mode = "markers", marker = list(size = 10, opacity = 0.7)) %>%
      layout(hovermode = "closest", title = "Interactive UMAP Plot") %>%
      layout(coloraxis = list(colorscale = "Viridis"), hoverinfo = "text")
  })
}

shinyApp(ui = ui, server = server)




















#converting .tif images to .png images so they will be read by the browser and displayed.

install.packages("magick")
library(magick)

tif_path = "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/nucleiSegmented/individuals/tif/"
png_path = "/Volumes/fsmresfiles/Basic_Sciences/CDB/CDB_Collaborations/Arispe_Goyal/MadelineMelzer/DATA/nucleiDysmorphia/data/age/nucleiSegmented/individuals/png/"


# List all TIFF files in the directory
tif_files <- list.files(tif_path, pattern = "\\.tif$", full.names = TRUE)

# Loop over each TIFF file and convert to PNG
for (tif_file in tif_files) {
  # Read the TIFF image
  img <- image_read(tif_file)

  # Define the output filename (replace .tif with .png)
  png_file <- file.path(png_path, paste0(basename(tools::file_path_sans_ext(tif_file)), ".png"))

  # Convert and save the image as PNG
  img %>% image_convert(format = "png") %>% image_write(png_file)
  # Optional: print a message to show progress
  cat(paste("Converted", tif_file, "to", png_file, "\n"))

}
