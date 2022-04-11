
library(shinyWidgets)
library(circlize)
library(dplyr)
library(factoextra)
library(cluster)
library(ggfortify)
library(ggrepel)
library(ComplexHeatmap)
library(shinythemes)
library(bslib)

#Load data
#link <- getURL("https://github.com/ouyaqing/Chemokine/blob/main/Transpose_data_V3.rds")
#link <- getURL("https://github.com/ouyaqing/Chemokine/blob/main/Brain_Chemo_Mus_analytics_ready.tsv")
#c.data <- read_delim(link, "\t", escape_double = FALSE, trim_ws = TRUE)
c.data <- readRDS("Transpose_data_V3.rds")

Category <- t(c.data[which(rownames(c.data) == 'Category'),])
Colour <- t(c.data[which(rownames(c.data) == 'Colour'),])

#Define UI

ui <- fluidPage(
tabsetPanel(
    id = "parameter",
    tabPanel(
      title = 'Parameters',
      pageWithSidebar(
        titlePanel("Interactive Chemokine-related Expression Analysis"),
        sidebarPanel(
          #select species
          radioButtons("species", label = "Choose species", 
                       choices = list("Human" = "Human","Mouse" = "Mouse"),
                       selected = "Human"),
          #select tissue
          radioButtons(inputId = "tissue",
                       label = "Choose tissues",
                       choices = list("Brain" = "Brain", "Lung" = "Lung", "Spleen" = "Spleen","Lymph-node" = "Lymph-node",'All tissues'=0),
                       selected = 0),
          #select analysing gene categories
          pickerInput(inputId = "genes",
                      label = "Choose genes",
                      choices = list("Chemokine receptors & ligands" = "Chemo","GAGs" = "GAGs","MMPs" ="MMPs"),
                      selected = "Chemo",
                      options = list('actions-box' = TRUE),
                      multiple = TRUE),
          #select cluster number
          sliderInput('clusters', 'Cluster count',value = 3,min = 1, max = 6),
          actionButton("pca", "Show the PCA plot",width = 200),
          actionButton("heatmap", "Show the correlation plot",width = 200)
        ),
        mainPanel(
          plotOutput('kmeans'),
          h4(textOutput("assay_var"),style = "padding-left:40px;padding-top:20px")
          
        )
      )
    ),
    tabPanel(
      title = "PCA",
      fluidRow(
        column(2,h4(textOutput("para_var1_1"))),
        column(2,h4(textOutput("para_var1_2"))),
        column(4,h4(textOutput("para_var1_3"))),
        column(7,plotOutput('plot_pca',width = 900,height = 800)),
        column(4,offset = 1,h4("Gene Lists:")),
        column(4,offset = 1,textOutput("gene_var1")),
        column(4,offset = 1,actionButton("return1", "Return",width = 160),style = "padding-top:20px"),
        

    )),
    tabPanel(
      title = "Correlation Matrix",
      fluidRow(
        column(2,h4(textOutput("para_var2_1"))),
        column(2,h4(textOutput("para_var2_2"))),
        column(4,h4(textOutput("para_var2_3"))),
        column(7,plotOutput('plot_heatmap',width = 700,height = 700)),
        column(4,offset = 1, h4("Gene Lists:")),
        column(4,offset = 1,textOutput("gene_var2")),
        column(4,offset = 1,actionButton("return2", "Return",width = 160),style = "padding-top:20px"),
        column(4,offset = 1,downloadButton("downloadData", "Download numeric table"),style = "padding-top:20px")
    
      )
    )
  )
)


#Define server funtion
server <- function(input, output,session) {
  #subset data
  # I.subset based on selected species & tissues & genes
    
    Final_data <- reactive({
      species_data <- subset(c.data,as.character(c.data$Species) %in% input$species)
      
      tissue_list <- c()
      row_names <- c()
      if (input$tissue != 0) {
        for (tissue in input$tissue){
          row_names <- rownames(species_data[species_data[,tissue]=='TRUE',])
          tissue_list <- c(tissue_list,row_names)
        }} else {
          tissue_list <- rownames(species_data)}
      tissue_list <- unique(tissue_list)
      tissue_data <- subset(species_data, rownames(species_data) %in% tissue_list)
      tissue_data <- as.data.frame(t(tissue_data))
      tissue_data <- cbind(tissue_data,Category,Colour)
      
      gene_data <- subset(tissue_data,as.character(tissue_data$Category) %in% input$genes)

      gene_data <- select(gene_data, -Category)
      
  #II. QC
      New_data <- gene_data
      # a. filter rows(genes) coverage > 70% assays (NA < 30)
      project_coverage <- round(rowSums(apply(is.na(New_data),2,as.numeric))/(ncol(New_data)-1)*100,1)
      New_data$Project_Cov <- project_coverage
      New_data_70 <- filter(New_data,Project_Cov < 30)
      New_data_70 <- select(New_data_70,-Project_Cov)
      New_data_90 <- as.data.frame(t(New_data_70))
      
      # b. filter column(assays) coverage > 90% genes(NA < 10)
      gene_coverage <- round(rowSums(apply(is.na(New_data_90),2,as.numeric))/(ncol(New_data_90)-1)*100,1)
      New_data_90$Gene_Cov <- gene_coverage
      New_data_90 <- filter(New_data_90,Gene_Cov < 10)
      New_data_90 <- select(New_data_90,-Gene_Cov)
      New_data_60 <- as.data.frame(t(New_data_90))
      
      # c. filter rows(genes) coverage > 60% assays (NA < 60)
      project_coverage2 <- round(rowSums(apply(is.na(New_data_60),2,as.numeric))/(ncol(New_data_60)-1)*100,1)
      New_data_60$Project_Cov <- project_coverage2
      New_data_60 <- filter(New_data_60,Project_Cov < 30)
      New_data_60 <- select(New_data_60,-Project_Cov)
      New_data_clean <- as.data.frame(t(New_data_60))
      
      # d. remove columns(assays) who have NA
      New_data_clean <- na.omit(New_data_clean)
      
      # e. remove coverage percentage column and row
      New_data_ready <- as.data.frame(t(New_data_clean)) 
      return(New_data_ready)
      
    })
  
    
  #load cluster number 
  cluster_num <- reactive({input$clusters})
  
  #plot the PCA
  output$kmeans <- renderPlot(fviz_nbclust(Final_data()[,1:(ncol(Final_data())-1)] %>% mutate_all(as.numeric), kmeans, method='silhouette'),width = 600,height = 400)
  
  #colour genes based on its annotation
  output$plot_pca <- renderPlot({
    autoplot(clara(Final_data()[,1:(ncol(Final_data())-1)]%>% mutate_all(as.numeric), cluster_num()),
             frame = T,
             size = 4,
             shape = 2, 
             fill.colour = Final_data()$Colour) + 
      geom_text_repel(aes(label = rownames(Final_data())), 
                      family="Courier",fontface ="bold", 
                      colour = Final_data()$Colour) + 
      theme(text = element_text(family="Times",face = "bold"))
      })
  
  #plot heatmap
  col_fun = colorRamp2(c(-1,-0.5, 0, 0.5,1), c('orange1',"gold", "white", 'tomato1',"orangered4"))
  heatmap_data <- reactive({
    df_transpose <- t((Final_data()[,1:(ncol(Final_data())-1)]%>% mutate_all(as.numeric)))
    Trans_data <- matrix(as.numeric(unlist(df_transpose)),nrow=nrow(df_transpose))
    colnames(Trans_data) <- rownames(Final_data())
    res <- cor(Trans_data, method = 'pearson') #correction matrix
    return(res)
    })
  
  output$plot_heatmap <- renderPlot({
      Heatmap(heatmap_data(), name = "Pearson",col=col_fun,row_names_gp = gpar(fontfamily='Times'),column_names_gp = gpar(fontfamily='Times'))
    })
  
  #buttons
  observeEvent(input$pca,{updateTabsetPanel(session,"parameter",selected = 'PCA')}) 
  observeEvent(input$heatmap,{updateTabsetPanel(session,"parameter",selected = 'Correlation Matrix')})
  observeEvent(input$return1,{updateTabsetPanel(session,"parameter",selected = 'Parameters')}) 
  observeEvent(input$return2,{updateTabsetPanel(session,"parameter",selected = 'Parameters')}) 

  #test output
  output$assay_var <- renderText({ 
    paste("Dataset: ",ncol(Final_data())-1, "assays are collected")
  })
  
  output$gene_var1 <- renderText({ 
    paste(rownames(Final_data()),collapse="; ")
  })
  
  output$gene_var2 <- renderText({ 
      paste(rownames(Final_data()),collapse="; ")
    })
  
  
  output$para_var1_1<- renderText({paste0("Species: ", input$species)})
  output$para_var1_2<- renderText({
    if (input$tissue != 0) {
      paste0("Tissue: ", input$tissue)}
    else{paste0("Tissue: All")
      }})
  output$para_var1_3<- renderText({paste0("Genes: ", list(input$genes))})
  output$para_var2_1<- renderText({paste0("Species: ", input$species)})
  output$para_var2_2<- renderText({
    if (input$tissue != 0) {
      paste0("Tissue: ", input$tissue)}
    else{paste0("Tissue: All")
    }})
  output$para_var2_3<- renderText({paste0("Genes: ", input$genes)})
  

  #download excel
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(Sys.Date(),"Correlation-Matrix",input$species,".csv", sep = "_")
    },
    content = function(file) {
      write.csv(heatmap_data(), file, row.names = TRUE)
    }
  )
}

shinyApp(ui = ui, server = server)
