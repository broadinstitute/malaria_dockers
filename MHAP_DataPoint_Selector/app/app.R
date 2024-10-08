library(shiny)
library(seqinr)

# Define UI
ui <- fluidPage(
  titlePanel("Data selector for MHap Analysis"),
  
  HTML("<div>
          <p>This application facilitates the selection of data points for MHap Analysis within the Terra platform. Its function is to compile a list of ASVs for inclusion in the analysis report and to define the time frame for analysis.</p>
          
          <p>To obtain the metadata file, follow these steps:</p>
          
          <ol>
              <li><strong>Upload FASTA Files:</strong> Begin by uploading the reference FASTA files using the Fasta Sequence Name Extractor. This action triggers the dynamic generation of checkboxes corresponding to each ASV name extracted from the uploaded files.</li>
              
              <li><strong>Upload Metadata File:</strong> Next, upload the metadata file using the Metadata File Reader. This step, similar to the previous one, dynamically loads the data in the metadafile. From the drop down menu, select the columns containing the time stamps. This will generate the checkboxes representing each time period relevant to the analysis.</li>
          </ol>
          
          <p>Once the checkboxes for genes and time stamps are populated, select the desired checkboxes and proceed to click the 'Download file' button. The generated file can then be deposited into the Terra workspace.</p>
      </div>"),
  tags$div(style="color: red;", HTML("<h4>As a minimum, one target ASV and one time stamp must be selected. The 'Download file' option will appear after this condition has been met.</h4>")),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("fasta_files", "Upload FASTA files", multiple = TRUE, accept = c(".fasta")),
      actionButton("select_all_seqs", "Select All"),
      actionButton("deselect_all_seqs", "Deselect All"),
      uiOutput("sequence_checkboxes"),
      br(),
      fileInput("csv_file", "Upload Metadata CSV file", accept = c(".csv")),
      selectInput("csv_column", "Choose a column:", choices = NULL),
      actionButton("select_all_entries", "Select All"),
      actionButton("deselect_all_entries", "Deselect All"),
      uiOutput("csv_checkboxes"),
      
      # Download button moved to below the sidebarPanel
      br(),
      uiOutput("download_button_ui")
    ),
    
    mainPanel(
      fluidRow(
        column(12, align = "center",
               tableOutput("selected_sequences"),
               tableOutput("selected_times")
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Function to extract sequence names from FASTA files
  extract_sequence_names <- function(files) {
    sequence_names <- NULL
    
    for (i in seq_along(files$datapath)) {
      fasta_data <- read.fasta(files$datapath[i])
      sequence_names <- c(sequence_names, names(fasta_data))
    }
    
    return(unique(sequence_names))
  }
  
  observeEvent(input$select_all_seqs, {
    updateCheckboxGroupInput(session, "selected_seqs", choices = unique(unlist(lapply(input$fasta_files$datapath, function(file_path) {
      fasta_data <- read.fasta(file_path)
      names(fasta_data)
    }))), selected = unique(unlist(lapply(input$fasta_files$datapath, function(file_path) {
      fasta_data <- read.fasta(file_path)
      names(fasta_data)
    }))))
  })
  
  observeEvent(input$deselect_all_seqs, {
    updateCheckboxGroupInput(session, "selected_seqs", choices = unique(unlist(lapply(input$fasta_files$datapath, function(file_path) {
      fasta_data <- read.fasta(file_path)
      names(fasta_data)
    }))), selected = character(0))
  })
  
  observeEvent(input$select_all_entries, {
    updateCheckboxGroupInput(session, "selected_csv_entries", choices = unique(unlist(lapply(input$csv_file$datapath, function(file_path) {
      df <- read.csv(file_path)
      unique(df[[input$csv_column]])
    }))), selected = unique(unlist(lapply(input$csv_file$datapath, function(file_path) {
      df <- read.csv(file_path)
      unique(df[[input$csv_column]])
    }))))
  })
  
  observeEvent(input$deselect_all_entries, {
    updateCheckboxGroupInput(session, "selected_csv_entries", choices = unique(unlist(lapply(input$csv_file$datapath, function(file_path) {
      df <- read.csv(file_path)
      unique(df[[input$csv_column]])
    }))), selected = character(0))
  })
  
  # Update checkboxGroupInput with sequence names
  output$sequence_checkboxes <- renderUI({
    req(input$fasta_files)
    
    sequence_names <- unlist(lapply(input$fasta_files$datapath, function(file_path) {
      fasta_data <- read.fasta(file_path)
      names(fasta_data)
    }))
    
    checkboxGroupInput("selected_seqs", "Select ASVs:", choices = unique(sequence_names))
  })
  
  # Read CSV file and update selectInput choices
  observe({
    req(input$csv_file)
    df <- read.csv(input$csv_file$datapath)
    updateSelectInput(session, "csv_column", choices = names(df))
  })
  
  # Update checkboxGroupInput with entries from selected CSV column
  output$csv_checkboxes <- renderUI({
    req(input$csv_file)
    req(input$csv_column)
    
    df <- read.csv(input$csv_file$datapath)
    choices <- unique(df[[input$csv_column]])
    
    checkboxGroupInput("selected_csv_entries", "Selected Time Stamps:", choices = choices)
  })
  
  # Render selected sequence names
  output$selected_sequences <- renderTable({
    req(input$selected_seqs)
    data.frame(Selected_ASVs = input$selected_seqs)
  })
  
  # Render selected times
  output$selected_times <- renderTable({
    req(input$selected_csv_entries)
    data.frame(Selected_Time_Stamps = input$selected_csv_entries)
  })
  
  # Render download button UI conditionally
  output$download_button_ui <- renderUI({
    if (!is.null(input$selected_seqs) && length(input$selected_seqs) > 0 &&
        !is.null(input$selected_csv_entries) && length(input$selected_csv_entries) > 0) {
      tags$div(style="width:200px;", downloadButton("download_selected", "Download file", style = "font-size: 20px; background-color: #4CAF50;"))
    }
  })
  
  # Generate downloadable CSV file
  output$download_selected <- downloadHandler(
    filename = function() {
      "mhap_data_points.csv"
    },
    content = function(file) {
      selected_checkboxes <- isolate(c("asvs__", input$selected_seqs, "timestamps__", input$selected_csv_entries))
      write.csv(data.frame(Selected_Checkboxes = selected_checkboxes), file, row.names = FALSE, quote = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
